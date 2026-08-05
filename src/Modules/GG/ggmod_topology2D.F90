!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides functionality to determine the 2D topological 
! mesh that is defined by the vessel and magnetic field. Although the 
! topological mesh is also a grid, it does not inherit from the GridUDT
! structure, since this structure typically contains too much or even
! irrelevant data. This does imply some duplication of some basic 
! routines, but that shouldn't be a big issue. The topological mesh can
! be used to generate an unstructured grid later on based on the 
! cells and flux tubes determined here. 

module ggmod_topology2D

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_dynamicarrays
    use mod_contour2D
    use mod_polygon
    use mod_graph, only: GraphUDT
    use mod_sort
    use mod_definitions
    use mod_linearsolverinterface, only: SolveDenseLinearSystemDI
    use goatmod_types, only : magneticFieldUDT, VesselUDT, &
        ConstructVesselPolygonSet
    use goatmod_userinput, only : TopomeshOptionsUDT
    use Interpolant1D
    use mod_streamlinetracing2D
    use PolygonLevelsetFunction2D
    use mod_sparseinterface
    use mod_constants, only: posinfval_R8
    use mod_utility, only: wall_time
    implicit none
    private 
    public :: TopomeshUDT, ConstructTopologicalMesh, TraceExtrema2D, &
        TraceTangencyPoints2D, ReadTopologicalMesh, WriteTopologicalMesh, &
        IdentifyTopologicalMeshType

    ! Module parameters
    real(R8), parameter, private        :: tprelfieldtol = 1e-14 ! relative field tolerance under which extrema are removed
    real(R8), parameter, private        :: disttol = 1e-12 ! distance tolerance
    real(R8), parameter, private        :: distfrac = 1e-3 ! distance fraction when removing edge vertices
    real(R8), private                   :: ts, te ! timing parameters
    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Topological vertex 
    type :: TopomeshVertUDT 

        ! Description
        !============
        ! Vertex structure for the topological mesh. Contains different 
        ! fields to identify the type of vertex and to navigate the mesh
        ! more easily. 

        integer(I8)                     :: ntot ! total number of vertices
        integer(I8), allocatable        :: ID(:), face(:), cell(:), &
            flags(:), type(:), faceP(:, :), cellP(:, :), fsID(:)
        real(R8), allocatable           :: x(:), y(:), fval(:)
        logical, allocatable            :: BV(:)

    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshVertex

        ! Getters
        procedure :: GetFace        => GetTMVertFace 
        procedure :: GetFaceNeig    => GetTMVertFaceNeig

    end type 

    ! Topological face
    type :: TopomeshFaceUDT

        ! Description
        !============
        ! Face structure for the topological mesh. Contains different 
        ! fields to identify the type of vertex and to navigate the mesh
        ! more easily. Note that contrary to classical grid faces, these
        ! faces are not a straight line between vessel vertices, but 
        ! contain a (simple) polygon that describes the spatial layout
        ! of the face. 

        integer(I8)                     :: ntot, nfc ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:, :), cell(:), &
            label(:), fsID(:), type(:), cellP(:, :)
        type(RealDynamicArrayUDT), allocatable  :: x(:), y(:)
        type(PolygonUDT), allocatable   :: pol(:)
        logical, allocatable            :: BF(:)

    contains 
    
        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshFace

        ! Getters
        procedure :: GetCell        => GetTMFaceCell

    end type

    ! Topological cell
    type :: TopomeshCellUDT

        ! Description
        !============
        ! Topological mesh cell

        integer(I8)                     :: ntot, ntotv, ntotf ! total number of vertices
        integer(I8), allocatable        :: ID(:), vert(:), vertP(:, :), &
            face(:), faceP(:, :), flags(:), tube(:)
    contains 

        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshCell

        ! Getters
        procedure :: GetVert        => GetTMCellVert
        procedure :: GetFace        => GetTMCellFace 

    end type

    ! Topological tubes
    type :: TopomeshTubeUDT

        ! Description
        !============
        ! Topological mesh tube object. This contains cells that are in 
        ! the same flux tube. The faces that are in this flux tube are 
        ! topological mesh faces that have the tube cells as boundaries.
        ! Note that these tubes can only be constructed as a 
        ! post-processing step after topological mesh construction. The
        ! cell faces and cells present in the 'cell' and 'face' structure
        ! are sorted after construction. Additionally, we store the 
        ! aligned faces and vertices from both sides in separate arrays,
        ! in this case in bndface1 (with pointer bndface1P), bndface2, &
        ! bndvert1, bndvert2. It should be noted that these faces and 
        ! vertices are not sorted in any particular direction (this is 
        ! non-trivial, since tube boundaries may be branching polygons)

        ! Note: a directional graph is now available after tube 
        ! construction that points in downward psi direction. This graph
        ! is used when determining which tubes to delete for example. 
        ! Other uses may be added in the future. Note that this graph is 
        ! only up to date after calling the topomesh interconnection    
        ! routine! 

        integer(I8)                                 :: ncell, nface, ntot
        integer(I8), allocatable, dimension(:)      :: cell, face,  &
            bndf1, bndf2, bndv1, bndv2, ftneig1, ftneig2, hfside
        integer(I8), allocatable, dimension(:, :)   :: cellP, faceP, &
            bndf1P, bndf2P, bndv1P, bndv2P, ftneig1P, ftneig2P
        logical, allocatable, dimension(:)          :: isclosed 
        type(GraphUDT)                              :: graph 

    contains
    
        ! Initializer
        procedure :: Initialize     => InitializeTopologicalMeshTube

        ! Deallocator
        procedure :: Deallocate     => DeallocateTopologicalMeshTube

        ! Getter
        procedure :: GetFace        => GetTMTubeFace
        procedure :: GetCell        => GetTMTubeCell
        procedure :: GetBndFace     => GetTMTubeBndFace
        procedure :: GetBndVert     => GetTMTubeBndVert
        procedure :: GetNeig        => GetTMTubeNeig
        procedure :: GetHighFluxNeig    => GetTMTubeHighFluxNeig
        procedure :: GetLowFluxNeig     => GetTMTubeLowFluxNeig
        procedure :: GetHighFluxBndFace     => GetTMTubeHighFluxBndFace
        procedure :: GetLowFluxBndFace      => GetTMTubeLowFluxBndFace
        procedure :: GetHighFluxBndVert     => GetTMTubeHighFluxBndVert
        procedure :: GetLowFluxBndVert      => GetTMTubeLowFluxBndVert

    end type 

    ! General topological mesh type
    type :: TopomeshUDT 

        type(TopomeshVertUDT)   :: vert 
        type(TopomeshFaceUDT)   :: face 
        type(TopomeshCellUDT)   :: cell 
        type(TopomeshTubeUDT)   :: tube 
        integer(I8)             :: nfs 
        type(IntegerDynamicArrayUDT)    :: fsID
        type(RealDynamicArrayUDT)       :: fsfval

    contains 

        ! Initialize
        procedure :: Initialize         => InitializeTopologicalMesh 

        ! Getters
        procedure :: GetInternalFaceIDs, GetBoundaryFaceIDs, &
            GetTargetFaceIDs, GetSeparatrixFaceIDs, GetVesselFaceIDs, &
            GetCoreFaceIDs, GetCoreCellIDs, GetWideGridCellIDs, &
            GetSeparatrixFluxSurfaceIDs, GetStrikePointIDs, GetXPOintIDs, &
            GetOPointIDs, GetPrimaryXPointIDs, GetStrikePointXPointIDs, &
            GetClosedContourTangencyPointIDs, GetLastFluxSurfaceFaceIDs
    end type 

    ! Topomesh adaptor 
    type :: TopomeshAdaptorUDT

        ! Description
        !============
        ! This object should be used to do more complicated topological
        ! mesh adaptations that require additional data to be passed
        ! between different structures. In particular, this object 
        ! stores rather expensive to compute data, such as psi 
        ! distributions on faces, and updates it accordingly when parts
        ! are removed to keep the computational cost (more) reasonable
        ! of the topomesh adaptations. Other quantities that are 
        ! modified and used in different adaptation routines can also be 
        ! stored here to simplify data handling. 

        ! Note: in practice, this type is primarily used for tube 
        ! merging. 

        real(R8)                                :: dpsimin, lradmin, &
            lpolsepmin 
        real(R8), allocatable, dimension(:)     :: tubedpsi, tubelrad
        integer(I8), allocatable, dimension(:)  :: illegalfsIDs, splitfsIDs
        type(RealDynamicArrayUDT), allocatable, dimension(:)    :: &
            facepsi, facedlcrad, facepsi_mono
        logical                                 :: allowsepmerge, &
            allowcoremerge, allowpfmerge, allowsepsimplify
        character(:), allocatable               :: mergetubemeth

        ! Store locally for convenience
        class(ContourtracerUDT), allocatable    :: fieldtracer
        class(StreamlinetracerUDT), allocatable :: streamlinetracer
        type(MagneticFieldUDT)                  :: magneticField 
        type(VesselUDT)                         :: vessel

    contains 

        ! Initialization routine
        procedure   :: Initialize               => InitializeTopomeshAdaptor

        ! Updaters
        procedure   :: RemoveFaceDataLogical    => RemoveFaceDataLogicalTA
        generic     :: RemoveFaceData           => RemoveFaceDataLogical
        procedure   :: AddFaceData              => AddFaceDataTA

        ! Topomesh tube merging routines
        procedure   :: EvaluateTMTubesMergeCriterion    => EvaluateTMTubesMergeCriterionTA
        procedure   :: MergeTMTubesDriver   => MergeTMTubesDriverTA ! main 
        procedure   :: MergeTMTubesSimple   => MergeTMTubesSimpleTA ! simple merger
        procedure   :: MergeTMTubesComplex  => MergeTMTubesComplexTA ! complex merger
        procedure   :: MergeTMTubes         => MergeTMTubesTA ! general merge selector
        procedure   :: MergeTMTubesOO       => MergeTMTubesOOTA ! open-open merge
        procedure   :: MergeTMTubesOC       => MergeTMTubesOCTA ! open-closed merge
        procedure   :: MergeTMTubesCC       => MergeTMTubesCCTA ! closed-closed merge
        procedure   :: MergeTMTubesS        => MergeTMTubesSTA ! separatrix (or other branching polygon) merge
        procedure   :: CollapseTMSeparatrixTube     => CollapseTMSeparatrixTubeTA ! separatrix tube collapser
        procedure   :: SplitTMTubes         => SplitTMTubesTA ! tube splitter
        procedure   :: GetMergeTubePairs    => GetMergeTubePairsWrapper ! merge tube pair getter
        procedure   :: GetMergeTubePairsExtensive   => GetMergeTubePairsExtensiveTA
        procedure   :: GetMergeTubePairsLocal       => GetMergeTubePairsLocalTA
        procedure   :: GetMergeTubeCycles   => GetTopomeshTubeCyclesTA ! merge cyclic tube getter
        procedure   :: CleanSplitTMTubes    => CleanSplitTMTubesTA

        ! Separatrix simplification
        procedure :: SimplifySeparatrices   => SimplifyTMSeparatricesTA 

        ! Auxiliary routines for tube merging
        procedure   :: IsTubePairMergeable  => IsTubePairMergeableTA
        procedure   :: IsTubeSplittable     => IsTubeSplittableTA

    end type

    contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                  TOPOLOGICAL MESH CONSTRUCTION                   !
    !------------------------------------------------------------------!

    ! Main constructor
    subroutine ConstructTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine generates the topological mesh based on the given magnetic
        ! field and vessel geometry. 

        ! Algorithm
        !==========
        ! 1) Compute points of interest (minima, maxima, saddle points, tangency
        ! points) 
        ! 2) Create a reduced topological mesh of the extrema to generate the cuts
        ! (should contain vertices and faces, not cells) - this forms the initial
        ! basis of the final topological mesh
        ! 3) Find the vessel boundary, introduce the tangency points and compute
        ! additional contours starting from those tangency points.
        ! 4) Compute all intersections between the current faces, introduce these
        ! as new vertices (and keep track of type). Intersections in end points are
        ! assumed to be known already and not added twice
        ! 5) Construct the topological mesh cells

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(in)             :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Auxiliary
        type(VesselUDT)                         :: tmvessel

        ! Initialize
        !===========
        ! Copy the vessel structure since this may change due to 
        ! reconstruction of the vessel polygons etc while constructing
        ! the topological mesh (also small geometrical features might 
        ! change, but that should be negligible)
        tmvessel = vessel

        ! Construct basic mesh
        !=====================
        ! Check if we construct from scratch or load in a file
        if (options%readexistingTM) then 
            ! Read
            call ReadTopologicalMesh(topomesh, options%TMfilepath)

            if (options%readexistingtracers) then 
                ! Read from file
                call fieldtracer%Read(options%TMfieldtracerfilepath)
                call vesseltracer%Read(options%TMvesseltracerfilepath)
            else
                ! Update the field tracer based on topomesh data
                call UpdateTracersFromTopomesh(topomesh, fieldtracer, &
                    magneticField, tmvessel, options)
            end if 
        else
            ! Construct from scratch
            call ConstructBasicTopologicalMesh(tmvessel, magneticField, options, &
                topomesh, fieldtracer, vesseltracer, streamlinetracer)
        end if

        ! Apply adaptations
        !==================
        if (options%doadaptations) then 
            call ModifyTopologicalMesh(tmvessel, magneticField, options, &
                topomesh, fieldtracer, vesseltracer, streamlinetracer)
        end if 

        ! Clean
        !======
        call CleanTopologicalMesh(topomesh)

        ! Write
        !======
        call WriteTopologicalMesh(topomesh, 'topomesh')

    end subroutine

    ! Basic topological mesh constructor
    subroutine ConstructBasicTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine generates the topological mesh based on the given magnetic
        ! field and vessel geometry. 

        ! Algorithm
        !==========
        ! 1) Compute points of interest (minima, maxima, saddle points, tangency
        ! points) 
        ! 2) Create a reduced topological mesh of the extrema to generate the cuts
        ! (should contain vertices and faces, not cells) - this forms the initial
        ! basis of the final topological mesh
        ! 3) Find the vessel boundary, introduce the tangency points and compute
        ! additional contours starting from those tangency points.
        ! 4) Compute all intersections between the current faces, introduce these
        ! as new vertices (and keep track of type). Intersections in end points are
        ! assumed to be known already and not added twice
        ! 5) Construct the topological mesh cells

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(inout)          :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Auxiliary
        real(R8)                                :: dxfracmin, dyfracmin
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vf, xgv, ygv
        real(R8), parameter                     :: emptyR8(0)= 0
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp, &
            xe, ye, fe
        integer(I8)                             :: nv, ntp 
        integer(I8), allocatable, dimension(:)  :: typee, IDs
        integer(I8), parameter                  :: emptyI8(0) = 0
        type(VesselUDT)                         :: newvessel

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Initialize data
        dxfracmin = 1e-4_R8
        dyfracmin = 1e-4_R8

        ! Initialize topomesh 
        call topomesh%Initialize()

        ! Determine domain bounds based on vessel and magnetic field extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv]), maxval([xps, magneticField%interp%xgv])]
        yb = [minval([yps, magneticField%interp%ygv]), maxval([yps, magneticField%interp%ygv])]

        ! Construct a 2D structured grid for tracing (may be extended
        ! in the future for different grid types)
        nv = options%vresx*options%vresy
        allocate(xg(nv), yg(nv), xgv(options%vresx), ygv(options%vresy))
        call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx,  options%vresy, 0.0_R8, 0.0_R8)

        ! Compute extrema
        !================
        ! Only for mesh refinement later on
        call TraceExtrema2D(xe, ye, fe, typee, fieldtracer, magneticField, &
            options%fdonewton)

        ! Add extrema & related boundaries
        !=================================
        ! Tangency points & vessel geometry
        call AddTopologicalMeshTangencyPoints2(topomesh, vesseltracer, &
            vessel, magneticField, options)

        ! Write output
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftertp', .false.)
        end if 

        ! Construct refined grid based on tangency points and extrema
        ntp = count((topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        allocate(xtp(ntp), ytp(ntp), Ftp(ntp), IDs(ntp))
        xtp = pack(topomesh%vert%x, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        ytp = pack(topomesh%vert%y, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        Ftp = pack(topomesh%vert%fval, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        IDs = pack(topomesh%vert%ID, &
            (topomesh%vert%type == TMvertextp1ID) .or. (topomesh%vert%type == TMvertextp2ID))
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, [xtp, xe], [ytp, ye], 5, 5, dxfracmin, dyfracmin)  

        ! Evaluate magnetic field and vessel
        allocate(Vf(size(xg)))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
        call magneticField%interp%Evaluate(xtp, ytp, 0, 0, Ftp)

        ! Update the field tracer 
        fieldtracer = ConstructStructuredTracer(&
            reshape(Vf, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp, IDs, spread(.false., 1, size(xtp)), &
            fieldtracer%npmin, fieldtracer%npmax, fieldtracer%dl)
        
        !! Visualize by tracing contours
        !resc = 100
        !dv = (maxval(Vf) - minval(Vf))
        !cgv = [(k, k = 0, resc)]*(dv*0.90_R8)/real(resc, kind=R8) + minval(Vf) + dv*0.05
        !contours = fieldtracer%TraceContours(cgv)
        !allocate(pcontours(size(contours)))
        !!$omp parallel do default(shared)
        !do k = 1, size(contours)
        !    call pcontours(k)%Construct(contours(k)%x, contours(k)%y)
        !end do 
        !!$omp end parallel do
        !call tempps%Construct(pcontours)
        !call tempps%WriteData('mfcontours')

        ! Extrema
        call AddTopologicalMeshExtrema(topomesh, fieldtracer, &
            magneticField, options)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afterextrema')
        end if 

        ! Remove parts that do not lie inside the vessel
        newvessel = vessel
        call TrimTopologicalMesh(topomesh, magneticField, newvessel)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afterextrema_trimmed')
        end if 

        ! Process points
        !===============
        ! Set extrema with nearly identical values to be the same value
        do i = 1, topomesh%vert%ntot-1
            if (topomesh%vert%type(i) == TMvertexsaddleID) then 
                do j = i+1, topomesh%vert%ntot 
                    if (topomesh%vert%type(j) == TMvertexsaddleID) then 
                        if (abs(topomesh%vert%fval(i)-topomesh%vert%fval(j)) < options%ffieldtol) then 
                            topomesh%vert%fval(j) = topomesh%vert%fval(i)
                            topomesh%vert%fsID(j) = topomesh%vert%fsID(i)
                        end if 
                    end if 
                end do 
            end if 
        end do 

        ! Set type 2 tangency points with nearly identical values to be the same value
        do i = 1, topomesh%vert%ntot-1
            if (topomesh%vert%type(i) == TMvertextp2ID) then 
                do j = i+1, topomesh%vert%ntot 
                    if (topomesh%vert%type(j) == TMvertextp2ID) then 
                        if (abs(topomesh%vert%fval(i)-topomesh%vert%fval(j)) < options%ffieldtol) then 
                            topomesh%vert%fval(j) = topomesh%vert%fval(i)
                            topomesh%vert%fsID(j) = topomesh%vert%fsID(i)
                        end if 
                    end if 
                end do 
            end if 
        end do 

        ! Compute necessary contours
        !===========================
        ! Add contours and intersections
        call AddTopologicalMeshContours(topomesh, magneticField, newvessel, &
            fieldtracer, options)

        ! Vertex faces (preliminary, for garbage tangency point removal)
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftercontours')
        end if 

        ! Remove garbage tangency points
        call RemoveGarbageTangencyPoints(topomesh)

        ! Eliminate limiter-like configurations
        !======================================
        ! Eliminate
        call RemoveTopologicalMeshLimiterRegions(topomesh, magneticField, &
            fieldtracer, streamlinetracer, options)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afterlimiterremoval')
        end if 

        ! Add necessary data
        !===================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)        

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Write
        !======
        ! Overwrite vessel structure
        vessel = newvessel

        ! Write topological mesh
        call WriteTopologicalMesh(topomesh, 'topomesh_base')

        ! Write tracers
        call fieldtracer%Write('TMfieldtracer')
        call vesseltracer%Write('TMvesseltracer')

    end subroutine

    ! Main modifier
    subroutine ModifyTopologicalMesh(vessel, magneticField, options, &
        topomesh, fieldtracer, vesseltracer, streamlinetracer)

        ! Description
        !============
        ! This routine contains all possible modifications that may be
        ! applied to a basic topological mesh. Though it is not recommended, one might 
        ! try to apply these modifications to an already modified mesh 
        ! if desired. No guarantees on result though (most of these 
        ! operations are irreversible)

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT)                       :: topomesh
        type(VesselUDT), intent(inout)          :: vessel
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        class(ContourTracerUDT), allocatable, intent(inout)  :: vesseltracer, fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer

        ! Auxiliary
        integer(I8)                             :: nTMv
        integer(I8), allocatable, dimension(:)  :: IDTMv
        real(R8), allocatable, dimension(:)     :: xTMv, yTMv, FTMv
        logical, allocatable, dimension(:)      :: keepTMv
        type(TopomeshAdaptorUDT)                :: tmadaptor

        ! Compute additional contours
        !============================
        ! Core boundaries? 
        if (options%addcoreboundaries) then 
            call AddTopologicalMeshCoreBoundaries(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if 

        ! Write output
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftercoreboundaries', .false.)
        end if 

        ! 'PF' boundaries?
        if (options%addPFboundaries) then 
            call AddTopologicalMeshPFBoundaries(topomesh, magneticField, &
                vessel, fieldtracer, options)
        end if

        ! Write output
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afterpfboundaries', .false.)
        end if 


        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Compute additional interconnnection data
        !=========================================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_beforecells')
        end if 

        ! Remove regions
        !===============
        ! Merge tubes (before aligned vessel parts)?
        ! Use adaptor
        call tmadaptor%Initialize(topomesh, fieldtracer, streamlinetracer, &
            magneticField, vessel, options)

        call tmadaptor%MergeTMTubesDriver(topomesh, options)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftermerge1')
        end if 

        ! Insert aligned vessel parts?
        if (options%alignvesselparts) then 
            ! Update the field tracer for the adjusted topomesh (just 
            ! add pseudo saddle points, grid shouldn't be updated)

            ! Construct refined grid based on tangency points and extrema
            keepTMv = (topomesh%vert%type == TMvertextp1ID) &
                .or. (topomesh%vert%type == TMvertextp2ID)  &
                .or. (topomesh%vert%type == TMvertexmaxID)  & 
                .or. (topomesh%vert%type == TMvertexminID)  & 
                .or. (topomesh%vert%type == TMvertexsaddleID)  
            nTMv = count(keepTMv)
            allocate(xTMv(nTMv), yTMv(nTMv), FTMv(nTMv), IDTMv(nTMv))
            xTMv = pack(topomesh%vert%x, keepTMv)
            yTMv = pack(topomesh%vert%y, keepTMv)
            FTMv = pack(topomesh%vert%fval, keepTMv)
            IDTMv = pack(topomesh%vert%ID, keepTMv)

            ! Update the tracer
            fieldtracer%xs = xTMv
            fieldtracer%ys = yTMv 
            fieldtracer%vs = FTMv 
            fieldtracer%IDs = IDTMv
            fieldtracer%order = spread(0_I8, 1, nTMv)

            ! Insert aligned parts
            call InsertAlignedVesselParts(topomesh, fieldtracer, &
                magneticField, vessel, options)
        end if 

        ! Merge tubes (after aligned vessel parts)?
        call tmadaptor%Initialize(topomesh, fieldtracer, streamlinetracer, &
            magneticField, vessel, options)

        call tmadaptor%MergeTMTubesDriver(topomesh, options)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftermerge2')
        end if 

        ! Remove regions if desired
        call RemoveTopologicalMeshRegions(topomesh, vessel, options)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afterregion')
        end if 

    end subroutine

    ! Extrema addition
    subroutine AddTopologicalMeshExtrema(topomesh, fieldtracer, &
        magneticField, options)

        ! Description
        !============
        ! This routine is largely based on TraceExtrema2DBox, but adds the extrema
        ! and their connection lines to the topology mesh as vertices and faces.
        ! Note that the location of the extrema can be refined using the Newton
        ! solver, but that the connecting faces are simply taken as the gradient
        ! lines and are therefore less accurate. For a description on how the
        ! extrema are computed, see TraceExtrema2DBox. 

        ! IMPORTANT: it is assumed that the boundary faces are already
        ! present in the topological mesh! These are used to determine
        ! which 'radial' faces, resulting from the dfdx, dfdy contours, 
        ! should be retained in the topological mesh. Only those that 
        ! do not intersect the boundary are retained. This means that
        ! only radial faces that start and end in an extremum without
        ! intersecting the vessel boundary will be retained. If there
        ! are multiple connections between two extrema, then only the 
        ! connection with monotonously increasing field value is 
        ! retained. This is required to have a conforming topomesh for 
        ! grid generation later on (otherwise these radial lines would
        ! intersect with contour lines of saddle points etc)

        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: when computing the optimum with the Newton solver, it is assumed
        ! that no damping strategy is necessary since the initial point should be
        ! (very) close to the optimum. 

        ! Note 3: we pass the field tracer to reuse the grid already
        ! used before. However, we need to initialize the tracers for
        ! dfdx and dfdy again, since 1) the values change, 2) the 
        ! saddle points of the original field should not be present 
        ! (these are in fact often no saddle points of the derivatives!)

        ! Note 4: only segments of dfdx and dfdy that start or end in a
        ! topological mesh extremum (saddle, min, max) are added to the 
        ! topological mesh

        ! Note 5: checking monotonically increasing field values can be
        ! very tricky, since discretization effects and small discrepancies
        ! between tracers (e.g. how saddle points are dealt with) may 
        ! quickly lead to minor non-monotonous behavior (typically at 
        ! start and end points but sometimes also somewhere along the
        ! curve). Therefore, we first find all potential segments and 
        ! rank them later on. If multiple connections are found, the one
        ! that is 'most' monotonous is retained, the others removed. 
        ! Alternatively, one can try to introduce streamlines in between
        ! the different extrema (that would actually also be better for
        ! grid generation). 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        integer(I8)                             :: ngp, nfxc, nfyc, &
            nx, nvinit
        integer(I8), allocatable, dimension(:)  :: ts1, ts2, tt, typee, &
            vf1, vf2, teid, tsid, sortind, temps1, temps2, v1, v2, tcIDs
        real(R8)                                :: tempx, tempy, ltot
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy, xe, ye, fe, tsr1, tsr2, &
            tsrid, tempxint, tempyint, tfval, dtfval, frac, dx, dy, dl
        logical                                 :: conv, donewton, addsegment
        logical, allocatable, dimension(:)      :: keepind
        integer, allocatable, dimension(:, :)   :: adjacencymatrix
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(RealDynamicArrayUDT), allocatable  :: xfrda(:), yfrda(:), &
            fxpsrid(:), fypsrid(:) 
        type(IntegerDynamicArrayUDT), allocatable   :: fxpeid(:), &
            fxpsid(:), fypeid(:), fypsid(:), cIDs(:)
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT)                        :: tempc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:), allc(:)
        type(PolygonUDT)                        :: temppol
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer
        type(PolygonSetUDT)                     :: tempps

        ! Loop
        integer(I8)                             :: i, j, k, ec, cc

        ! Initialize
        !===========
        ! Check inputs
        donewton = options%fdonewton
        if (.not. donewton) then 
            print *, 'TraceExtrema2DBox: not applying Newton solver, ' // & 
                'locations of extrema may be slightly inaccurate'
        end if
        allocate(allc(0))

        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Initialize extrema counter
        ec = 0

        ! Store initial amount of vertices
        nvinit = topomesh%vert%ntot ! to update intersection IDs later on

        ! Trace contours
        !===============
        ! Evaluate derivatives
        allocate(f(ngp), fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, f)
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs, fxtracer%issaddle)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0), fxtracer%issaddle(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs, fytracer%issaddle)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0), fytracer%issaddle(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)
        

        ! Trace 
        fxc = fxtracer%TraceContours([0.0_R8])
        fyc = fytracer%TraceContours([0.0_R8])

        ! Clean
        call CleanContours(fxc)
        call CleanContours(fyc)

        ! Compute all intersections
        !==========================
        ! Convert all to polygons
        nfxc = size(fxc)
        nfyc = size(fyc)
        allocate(fxp(nfxc), fyp(nfyc))

        do i = 1, nfxc 
            call fxp(i)%Construct(fxc(i)%x, fxc(i)%y)
        end do 
        do i = 1, nfyc 
            call fyp(i)%Construct(fyc(i)%x, fyc(i)%y)
        end do 

        ! Write out data
        call tempps%Construct(fxp)
        call tempps%WriteData('extrema_fx_lines')
        call tempps%Construct(fyp)
        call tempps%WriteData('extrema_fy_lines')
        

        ! Initialize intersection trackers
        allocate(fxpeid(nfxc), fxpsid(nfxc), fypeid(nfyc), fypsid(nfyc), &
            fxpsrid(nfxc), fypsrid(nfyc))
        do i = 1, nfxc 
            fxpeid(i) = ConstructIntegerDynamicArray()
            fxpsid(i) = ConstructIntegerDynamicArray()
            fxpsrid(i) = ConstructRealDynamicArray()
        end do 
        do i = 1, nfyc 
            fypeid(i) = ConstructIntegerDynamicArray()
            fypsid(i) = ConstructIntegerDynamicArray() 
            fypsrid(i) = ConstructRealDynamicArray()
        end do 

        ! Compute intersections
        !$omp parallel do default(none) schedule(dynamic) collapse(2) ordered &
        !$omp shared(fxp, fyp, magneticField, fxpeid, fxpsid, fxpsrid, &
        !$omp fypeid, fypsid, fypsrid, xc, yc, fc, tc, ec, nfxc, nfyc, donewton) & 
        !$omp private(i, j, tx, ty, ts1, ts2, tsr1, tsr2, nx, tt, k, &
        !$omp tempx, tempy, conv, tf, tfxx, tfxy, tfyy, thiseig)
        do i = 1, nfxc
            ! Only compute intersections with other polygons
            do j = 1, nfyc
                ! Compute intersections with next polygon
                call PolygonIntersections(fxp(i), fyp(j), tx, ty, ts1, ts2, &
                    tsr1, tsr2)
                
                ! Check if found
                nx = size(tx)
            
                ! Add intersections
                allocate(tt(nx))
                if (donewton) then 
                    do k = 1, nx 
                        ! Refine
                        call TinyNewtonSolver(tempx, tempy, conv, &
                            tx(k), ty(k), magneticField)
                        
                        ! Check 
                        if (conv) then 
                            tx(k) = tempx
                            ty(k) = tempy
                        else
                            print *, 'TraceExtrema2DBox: newton solver did ' // & 
                                'not converge, taking original estimate'
                        end if
                    end do 
                end if 

                ! Compute value at location and second order derivatives
                allocate(tf(nx), tfxx(nx), tfxy(nx), tfyy(nx))
                call magneticField%interp%Evaluate(tx, ty, 0, 0, tf)
                call magneticField%interp%Evaluate(tx, ty, 2, 0, tfxx)
                call magneticField%interp%Evaluate(tx, ty, 1, 1, tfxy)
                call magneticField%interp%Evaluate(tx, ty, 0, 2, tfyy)

                ! Determine order
                do k = 1, nx

                    ! Compute value and check type of extremum
                    thiseig = ComputeEigenvaluesSymmetric2by2Matrix(&
                        tfxx(k), tfyy(k), tfxy(k))
                    if (any(thiseig == 0)) then 
                        print *, 'TraceExtrema2DBox: extrema with zero ' // & 
                            'eigenvalue detected, may be improperly identified'
                    end if 
                    if (all(thiseig > 0)) then 
                        ! Local minimum
                        tt(k) = -1;
                    elseif (all(thiseig < 0)) then 
                        ! Local maximum
                        tt(k) = 1;
                    else
                        ! Saddle point
                        tt(k) = 0;
                    end if

                end do

                ! Store IDs, connectivity, and values together so each
                ! extremum ID remains aligned with its coordinate-array index
                ! when contour pairs are processed by different threads.
                !$omp ordered
                do k = 1, nx
                    ec = ec + 1

                    ! Store intersection data
                    call fxpeid(i)%Append(ec)
                    call fxpsid(i)%Append(ts1(k))
                    call fxpsrid(i)%Append(tsr1(k))
                    call fypeid(j)%Append(ec)
                    call fypsid(j)%Append(ts2(k))
                    call fypsrid(j)%Append(tsr2(k))
                end do

                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)
                !$omp end ordered

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        !$omp end parallel do 
        
        ! Add extrema
        !============
        ! Extract
        xe = xc%Get()
        ye = yc%Get()
        fe = fc%Get()
        typee = tc%Get()

        ! Set vertex type
        do i = 1, size(typee)
            if (typee(i) == -1_I8) then 
                typee(i) = TMvertexminID 
            elseif (typee(i) == 0_I8) then 
                typee(i) = TMvertexsaddleID 
            else
                typee(i) = TMvertexmaxID 
            end if 
        end do 

        ! Loop to add
        do i = 1, ec 
            ! Update ID
            topomesh%nfs = topomesh%nfs+1

            ! Add vertex
            call AddTopologicalMeshVertex(topomesh, xe(i), ye(i), fe(i), &
                typee(i), topomesh%nfs)
            
        end do

        ! Reconstruct faces
        !==================
        ! Initialize contour counter
        cc = 0

        ! Loop over fxp
        do i = 1, nfxc 
            ! Unpack
            teid = fxpeid(i)%Get()
            tsid = fxpsid(i)%Get()
            tsrid = fxpsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)

            ! Check
            if (size(teid) == 0) then 
                ! No intersections whatsoever, so no saddle points at 
                ! start and end - continue without adding
                cycle
            end if 
            
            ! Extract the faces from this polygon
            call ExtractTopologicalFacesFromPolygon(fxp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Insert into topological mesh as contour (need to rebuild contour)
            tempc = fxc(i)
            do j = 1, size(xfrda)
                ! Check
                addsegment = .true. 
                if ((vf1(j) == 0) .or. (vf2(j) == 0)) then 
                    ! Skip
                    cycle
                end if 

                ! Check if we intersect a boundary face
                call temppol%Construct(xfrda(j)%Get(), yfrda(j)%Get())
                do k = 1, topomesh%face%ntot
                    if (topomesh%face%type(k) == TMfacebndID) then 
                        ! Check for intersections
                        call PolygonIntersections(topomesh%face%pol(k), &
                            temppol, tempxint, tempyint, temps1, temps2)
                        if (size(tempxint) /= 0) then 
                            ! Skip and exit
                            addsegment = .false. 
                            exit 
                        end if 
                    end if 
                end do 

                

                ! Check if we add
                if (addsegment) then 
                    tempc%x = xfrda(j)%Get()
                    tempc%y = yfrda(j)%Get()
                    tempc%startsaddle = vf1(j) + nvinit ! per definition non-zero, so this is fine
                    tempc%endsaddle = vf2(j) + nvinit
                    allc = [allc, tempc]
                    !call InsertTopologicalMeshContour(topomesh, magneticField, &
                    !    tempc, TMfaceradID, 0)
                end if 
            end do 

            ! Add to faces
            !do j = 1, size(vf1)
            !    call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
            !        xfrda(j), yfrda(j), TMfaceradID, 0)
            !end do 
        end do 

        ! Loop over fyp
        do i = 1, nfyc 
            ! Unpack
            teid = fypeid(i)%Get()
            tsid = fypsid(i)%Get()
            tsrid = fypsrid(i)%Get()
            allocate(sortind(size(tsrid)))
            call Sort(tsrid, ind=sortind)
            tsid = tsid(sortind)
            teid = teid(sortind)
            deallocate(sortind)

            ! Check
            if (size(teid) == 0) then 
                ! No intersections whatsoever, so no saddle points at 
                ! start and end - continue without adding
                cycle
            end if 
            
            ! Extract the faces from this polygon
            call ExtractTopologicalFacesFromPolygon(fyp(i), teid, &
                tsid, xe, ye, vf1, vf2, xfrda, yfrda)

            ! Insert into topological mesh as contour (need to rebuild contour)
            tempc = fyc(i)
            do j = 1, size(xfrda)
                ! Check
                addsegment = .true.
                if ((vf1(j) == 0) .or. (vf2(j) == 0)) then 
                    ! Skip
                    cycle
                end if 

                ! Check if we intersect a boundary face
                call temppol%Construct(xfrda(j)%Get(), yfrda(j)%Get())
                do k = 1, topomesh%face%ntot
                    if (topomesh%face%type(k) == TMfacebndID) then 
                        ! Check for intersections
                        call PolygonIntersections(topomesh%face%pol(k), &
                            temppol, tempxint, tempyint, temps1, temps2)
                        if (size(tempxint) /= 0) then 
                            ! Skip and exit
                            addsegment = .false. 
                            exit 
                        end if 
                    end if 
                end do 

                ! Check if we need to add 
                if (addsegment) then 
                    tempc%x = xfrda(j)%Get()
                    tempc%y = yfrda(j)%Get()
                    tempc%startsaddle = vf1(j) + nvinit ! per definition non-zero, so this is fine
                    tempc%endsaddle = vf2(j) + nvinit
                    allc = [allc, tempc]
                    !call InsertTopologicalMeshContour(topomesh, magneticField, &
                    !    tempc, TMfaceradID, 0)
                end if 
            end do 

        end do
        
        ! Check which extrema are interconnected
        v1 = [allc%startsaddle]
        v2 = [allc%endsaddle]
        allocate(adjacencymatrix(maxval([v1, v2]), maxval([v1, v2])))
        adjacencymatrix = 0
        cc = 0 ! connection counter
        do i = 1, size(v1)
            if (adjacencymatrix(v1(i), v2(i)) == 0 .and. adjacencymatrix(v2(i), v1(i)) == 0) then 
                cc = cc + 1
                adjacencymatrix(v1(i), v2(i)) = cc 
                adjacencymatrix(v2(i), v1(i)) = cc 
            end if 
        end do 

        ! Determine which parts connect each pair of extrema
        allocate(cIDs(cc), keepind(size(v1)))
        keepind = .false. 
        do i = 1, cc
            cIDs(i) = ConstructIntegerDynamicArray()
        end do
        do i = 1, size(v1)
            call cIDs(adjacencymatrix(v1(i), v2(i)))%Append(i)
        end do 

        ! For each pair, determine which segment to keep
        do i = 1, size(cIDs)
            ! Check
            if (cIDs(i)%Size() > 1) then 
                ! Get contour parts
                tcIDs = cIDs(i)%Get()

                ! Determine non-monotonous length fraction (smallest)
                allocate(frac(size(tcIDs)))
                frac = 0.0_R8
                do j = 1, size(tcIDs)
                    ! Get coordinates etc
                    tx = allc(tcIDs(j))%x
                    ty = allc(tcIDs(j))%y
                    dx = tx(2:) - tx(1:size(tx)-1)
                    dy = ty(2:) - ty(1:size(ty)-1)
                    dl = sqrt(dx**2 + dy**2)
                    ltot = sum(dl)
                    tfval = fieldtracer%Evaluate(tx, ty)
                    dtfval = tfval(2:) - tfval(1:size(tfval)-1)

                    ! Determine length fraction that is non-monotonous
                    if (tfval(1) > tfval(size(tfval))) then 
                        ! Decreasing, set edge lengths to zero where they
                        ! are aligned
                        where (dtfval >= 0.0_R8) dl = 0.0_R8
                    else
                        ! Increasing, set edge lengths to zero where they
                        ! are aligned
                        where (dtfval <= 0.0_R8) dl = 0.0_R8
                    end if 
                    frac(j) = sum(dl)/ltot
                    frac(j) = abs((maxval(tfval) - minval(tfval))/abs(tfval(1) - tfval(size(tfval))))
                end do 

                ! Keep the one where frac is minimal
                keepind(tcIDs(minloc(frac))) = .true.

                ! Housekeeping 
                deallocate(frac)
            else
                ! Mark this contour to be kept
                keepind(cIDs(i)%Get()) = .true. 
            end if
        end do

        ! Introduce contours in topomesh
        do i = 1, size(allc)
            if (keepind(i)) then 
                call InsertTopologicalMeshContour(topomesh, magneticField, &
                    allc(i), TMfaceradID, 0)
            end if 
        end do

    end subroutine 

    ! Tangency points and vessel addition
    subroutine AddTopologicalMeshTangencyPoints(topomesh, tpx, tpy, tpf, &
            tptype, bndtracer)

        ! Description
        !============
        ! This routine adds the tangency points and the vessel segments they
        ! construct to the topological mesh. To this end, the tangency points
        ! should already be traced using the TraceTangencyPoints2DBox routine. The
        ! vessel polygon is retraced (using resolution resx, resy) and the tangency
        ! points are then inserted at the nearest edge. Note that this may result
        ! in slight vessel geometry changes. 

        ! Algorithm
        !==========
        ! Actually it is quite easy: we add the tangency point as 'saddle' points,
        ! and make sure the value of those saddle points is exactly zero. When
        ! tracing the field line starting from one of the tangency points, the
        ! TraceContourLineSTructured2DPoint routine will return a single contour
        ! with all tangency points properly included. To make sure we trace each
        ! separate contour, we loop until all tangency points have been found. Easy
        ! peasy lemon squeazy

        ! Note: it is assumed that the tangency points were already 
        ! added to the boundary tracer beforehand! 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)   :: topomesh 
        class(ContourTracerUDT), intent(in) :: bndtracer
        real(R8), intent(in)                :: tpx(:), tpy(:), tpf(:)
        integer(I8), intent(in)             :: tptype(:)

        ! Auxiliary
        integer(I8)                         :: vcinit, ttype, ttp, &
            ftype, facevert(1:2), fsID
        logical, allocatable, dimension(:)  :: notfound 
        type(ContourUDT), allocatable       :: tc(:)
        type(RealDynamicArrayUDT)           :: tx, ty 

        ! Loop
        integer(I8)                 :: i 

        ! Initialize
        !===========
        ! Counters (relative to current state of topomesh)
        vcinit = topomesh%vert%ntot

        ! Add tangency points and construct mapping
        do i = 1, size(tpx)
            if (tptype(i) == 1) then 
                ttype = TMvertextp1ID 
            else
                ttype = TMvertextp2ID
            end if 
            call AddTopologicalMeshVertex(topomesh, tpx(i), tpy(i), &
                tpf(i), ttype, topomesh%nfs+1)
            topomesh%nfs = topomesh%nfs + 1
        end do

        ! Loop
        !=====
        allocate(notfound(size(tpx)))
        notfound = .true. 
        do while (any(notfound)) 
            ! Find first tangency point to trace from
            ttp = findloc(notfound, .true., 1)

            ! Sanity check
            if (ttp == 0) then 
                call gdErrorHandler('AddTopologicalMeshTangencyPoints: bug detected')
            end if 
            
            ! Trace the contour
            tc = bndtracer%TraceContours([tpx(ttp)], [tpy(ttp)])
            
            ! Add faces & set notfound
            do i = 1, size(tc)
                ! Check
                if ((tc(i)%startsaddle == 0) .or. (tc(i)%endsaddle == 0)) then 
                    call gdErrorHandler('AddTopologicalMeshTangencyPoints: ' // & 
                        'contour lines do not start and end in tangency ' // & 
                        'points for vessel, check input')
                end if 
                
                ! Make sure no segments are added that close upon themselves -
                ! shouldn't be possible at this stage, only in some 
                ! very exceptional cases
                if (tc(i)%startsaddle == tc(i)%endsaddle) then 
                    print *, 'AddTopologicalMeshTangencyPoints: segment detected that closes upon itself. Removing...'
                    cycle
                end if 

                ! Add the face
                ftype = TMfacebndID ! vessel boundary -> outer boundary
                facevert = [tc(i)%startsaddle, tc(i)%endsaddle];
                notfound(facevert) = .false.
                facevert = facevert + vcinit; ! adjust here to get the correct global vertex index!
                fsID = 0; ! no ID, because vessel boundary
                tx = ConstructRealDynamicArray(tc(i)%x)
                ty = ConstructRealDynamicArray(tc(i)%y)
                call AddTopologicalMeshFace(topomesh, facevert, tx, ty, ftype, fsID, 0.0_R8)
            end do 
        end do



    end subroutine

    ! 2D tangency points
    subroutine AddTopologicalMeshTangencyPoints2(topomesh, &
        boundarytracer, vessel, magneticField, options)

        ! Description
        !============
        ! This routine traces the tangency points on a prescribed boundary.
        ! We base ourselves solely on the
        ! discrete representation of the vessel geometry and compute the tangency
        ! points by looking at the field evaluated in the vertices of the vessel
        ! polygon(s). Here, minima and maxima are quite easily and rapidly found
        ! without having to compute any intersections between curves. To determine
        ! the type of tangency point (i.e. whether the curve bends inside or
        ! outside the domain), we check whether it's a minima or maxima along the
        ! vessel curve and whether the tangency point field value is found in the
        ! interior of the domain near the tangency point location. 

        ! The latter is a bit tricky here: since we don't explicitly know the
        ! magnetic field in the interior, we evaluate the derivative in the
        ! tangency point location and base ourselves on that. This, however, may be
        ! inaccurate since the tangency point location itself will not correspond
        ! with the location of the continuous optimum. Caution is therefore advised
        ! when the local magnetic field varies strongly locally (or when the mesh
        ! is too coarse). 


        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: it is assumed that the vessel polygon forms a closed (or a set of
        ! closed) surfaces. These polygons are oriented later on such that the face
        ! normal points inward the domain at all times. 


        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(in)     :: boundarytracer
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options
        type(VesselUDT), intent(in)             :: vessel

        ! Auxiliary
        integer(I8)                             :: flag, ntp
        integer(I8), allocatable                :: dvals(:), ddvals(:), &
            tv(:), extrlocind(:), tt(:), vf1(:), vf2(:), eID(:), sID(:), &
            sortind(:)
        real(R8)                                :: fdifftol 
        real(R8), allocatable, dimension(:)     :: val, dval, tx, ty, &
            tf, nxpe, nype, nxp, nyp, normprod, dFdx, dFdy
        logical                                 :: hasbeendeleted
        logical, allocatable                    :: extrloc(:), keepind(:), &
            doflip(:)
        type(RealDynamicArrayUDT), allocatable  :: xf(:), yf(:)
        type(ContourUDT), allocatable           :: bndcontours(:)
        type(PolygonUDT), allocatable           :: bndpol(:), realbndpol(:)
        type(PolygonSetUDT)                     :: bndps

        ! Loop
        integer(I8)                             :: i, j, k 

        ! Initialize
        !===========
        ! Initialize dynamic arrays

        ! Construct boundary polygonset
        !==============================
        ! Check if we just use the existing vessel polygon set or if
        ! we retrace 
        if (options%dotpvesselbased) then
            ! Unpack
            bndps = vessel%polygonset
            bndpol = bndps%polygons
        else
            ! Boundary contours
            bndcontours = boundarytracer%TraceContours([0.0_R8])
            
            ! Sanity check
            if (size(bndcontours) == 0) then 
                ! Throw error - no boundary polygon
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: could not trace ' // & 
                    'boundary contour, check input')
            end if 

            
            ! Construct polygons from boundary
            allocate(bndpol(size(bndcontours)), keepind(size(bndcontours)))
            keepind = .true. 
            do i = 1, size(bndcontours)
                ! Construct
                call bndpol(i)%Construct(bndcontours(i)%x, bndcontours(i)%y)

                ! Check for closedness, if not -> warning and remove
                if (.not. bndpol(i)%isclosed) then 
                    print *, 'AddTopologicalMeshTangencyPoints2: ' //&
                        'boundary polygon part is not closed and will not be ' // &
                        'added'
                    print *, 'part: ', i
                    keepind(i) = .false. 
                    !call gdErrorHandler('AddTopologicalMeshTangencyPoints2: boundary polygon ' // & 
                    !    'is not closed, not supported. Check input')
                end if 
            end do 

            ! Get only closed polygons
            allocate(realbndpol(count(keepind)))
            realbndpol = pack(bndpol, keepind)
            if (size(realbndpol) == 0) then 
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: no ' // & 
                    'closed boundary polygons found, check input')
            end if 
            bndpol = realbndpol

            ! Construct polygonset
            call bndps%Construct(bndpol)

            ! Orient
            call bndps%OrientNestedClosedPolygons(flag)

            ! Check if successful
            if (flag /= 0) then
                call gdErrorHandler('AddTopologicalMeshTangencyPoints2: could not orient ' // & 
                    'boundary polygons, check input')
            end if 
        end if

        ! Compute tangency points
        !========================
        do i = 1, size(bndps%polygons)
            ! Associate for ease
            associate(p         => bndps%polygons(i))

            ! Evaluate field at vertex locations
            allocate(val(size(p%vert)))
            call magneticField%interp%Evaluate(p%x(p%vert), p%y(p%vert), &
                0, 0, val) ! assumed start and end point the same
            val = [val, val(2)] ! extend to take next edge into account

            ! Compute precision
            fdifftol = (maxval(val) - minval(val))*tprelfieldtol ! why is that still in here?
            fdifftol = 0.0_R8 

            ! Take difference
            dval = val(2:size(val)) - val(1:size(val)-1)

            ! Check where this changes sign
            allocate(dvals(size(dval)))
            dvals = 0
            where (dval > fdifftol) dvals = 1
            where (dval < -fdifftol) dvals = -1
            ddvals = dvals(2:size(dvals)) - dvals(1:size(dvals)-1) 

            ! Find the location and value of extrema
            extrloc = [.false., ddvals /= 0]
            allocate(tv(count(extrloc)))
            tv = pack(p%vert, extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = pack(val(1:size(val)-1), extrloc)
            
            ! Check if we should exclude extremum pairs based on field value
            ! difference
            
            k = 1
            hasbeendeleted = .false.
            !do while (k < size(tf))
            !    ! Check difference
            !    if ((tv(k+1)-tv(k) == 1)) then ! (abs(tf(k+1)-tf(k)) < fdifftol) .or. 
            !        ! Remove values, such that subsequent ones can be
            !        ! checked too
            !        hasbeendeleted = .true.
            !        extrloc(tv(k+1)) = .false.
            !        tv = [tv(:k-1), tv(k+1:)]
            !        tf = [tf(:k-1), tf(k+1:)]
            !        k = 1
            !    else
            !        k = k + 1
            !    end if
            !end do 
            
            ! Check last 'edge' - shouldn't do this, since boundary
            ! polygons are closed and therefore first and last 
            ! vertex are the same
            !if (abs(tf(size(tf))-tf(1)) < fdifftol) then 
            !    extrloc(tv([1, size(tf)])) = .false.
            !    tv = pack(p%vert, extrloc)
            !    tf = val(tv)
            !    hasbeendeleted = .true.
            !end if
            
            ! Issue message
            if (hasbeendeleted) then 
                print *, 'AddTopologicalMeshTangencyPoints2: some ' // & 
                    'tangency points were deleted based on their field ' // & 
                    'values as they are very close together'
            end if

            ! Recompute points
            tv = pack(p%vert, extrloc)
            extrlocind = tv 
            extrlocind = pack([(k, k = 1, size(extrloc))], extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(extrlocind)

            ! Compute product between normal at vertex and magnetic field
            allocate(dFdx(size(tf)), dFdy(size(tf)))
            call magneticField%interp%Evaluate(tx, ty, 1, 0, dFdx)
            call magneticField%interp%Evaluate(tx, ty, 0, 1, dFdy)

            ! Get normals of edges
            nxpe = [p%nx, p%nx(1)] 
            nype = [p%ny, p%ny(1)]

            ! Check orientation - default polygon orientation is inward pointing,
            ! we need outward pointing
            doflip = [p%edges(:, 1) == p%vert(1:p%ne)] ! if equal, then the normal points inwards and we need to flip
            doflip = [doflip, doflip(1)]
            where (doflip)
                nxpe = -nxpe
                nype = -nype 
            end where

            ! Get normals at vertices
            nxp = 0.5*(nxpe(extrlocind-1) + nxpe(extrlocind))
            nyp = 0.5*(nype(extrlocind-1) + nype(extrlocind))

            ! Compute dot product
            normprod = dFdx*nxp + dFdy*nyp

            ! Determine type and add
            allocate(tt(size(tf)))
            tt = 0
            do j = 1, size(extrlocind)

                ! Determine type
                if ((val(extrlocind(j)-1) < tf(j)) .and. &
                    (val(extrlocind(j)+1) < tf(j))) then 

                    ! Local maximum
                    if (normprod(j) < 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = TMvertextp1ID
                    else
                        ! Curve bends inwards
                        tt(j) = TMvertextp2ID
                    end if
                else
                    ! Local minimum
                    if (normprod(j) > 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = TMvertextp1ID
                    else
                        ! Curve bends inwards
                        tt(j) = TMvertextp2ID
                    end if 
                end if
            end do 

            ! Add to topological mesh
            !========================
            ! Extract faces from polygon
            eID = [(k, k = 1, size(tv))]
            sID = tv 
            ntp = size(tv)
            allocate(sortind(size(sID)))
            call Sort(sID, ind=sortind)
            eID = eID(sortind)

            ! Check if the first vertex is a tangency point. If so, add 
            ! the first vertex as well
            if (sID(1) == 1) then 
                sID = [sID, bndpol(i)%ne]
                eID = [eID, eID(1)]
            end if 
            if (sID(size(sID)) == bndpol(i)%ne .and. (sID(1) /= 1)) then 
                sID = [1, sID]
                eID = [eID(size(eID)), eID]
            end if 
            call ExtractTopologicalFacesFromPolygon(bndpol(i), eID, sID, &
                tx, ty, vf1, vf2, xf, yf)  
            deallocate (sortind) 

            ! Update face vertices to correct indices for after adding
            ! vertices
            vf1 = vf1 + topomesh%vert%ntot 
            vf2 = vf2 + topomesh%vert%ntot 

            ! Add vertices
            do j = 1, ntp 
                call AddTopologicalMeshVertex(topomesh, tx(j), ty(j), &
                    tf(j), tt(j), topomesh%nfs+1)
                topomesh%nfs = topomesh%nfs + 1
            end do 

            ! Add faces
            do j = 1, size(xf)
                call AddTopologicalMeshFace(topomesh, [vf1(j), vf2(j)], &
                    xf(j), yf(j), TMfacebndID, 0_I8, 0.0_R8)
            end do 

            ! Housekeeping
            deallocate(val, dvals, dFdx, dFdy, tt, tv)
            end associate

        end do 
    end subroutine

    ! Necessary contours
    subroutine AddTopologicalMeshContours(topomesh, magneticField, vessel, &
        fieldtracer, options)

        ! Description
        !============
        ! This routine adds all necessary remaining contours of the magnetic field
        ! and computes their intersection(s) with the existing boundaries. It is
        ! assumed that all radial lines (cuts, vessel boundaries, etc) are already
        ! properly added to the topological mesh and that sufficient information is
        ! available for which points to compute the contours. The contours of the
        ! following points should be traced, starting from that point:
        !
        ! - tangency points with contour in the vessel (type 5)
        ! - saddle points (type 2)
        ! - split points on boundaries 
        !
        ! At this point, no contours for maxima/minima (types 1, 3) or for other
        ! tangency points (type 4) are added. This can be done in a later step for
        ! gridding purposes (i.e. if one does not want a region to be gridded near
        ! a maximum or tangency point).

        ! Algorithm
        !==========
        ! 1) Compute all additional contour parts starting from the points of
        ! interest. Make sure to parse all possible saddle point locations to the
        ! tracing routine.
        ! 2) Compute all intersections between all current boundaries and the new
        ! contours (contours themselves shouldn't intersect!)
        ! 3) Add the new vertices and faces to the topological mesh
        ! 4) Remove all points outside of the boundary and delete any associated
        ! faces, vertices, ...
        ! 5) Check for any faces that are closed or duplicate. If there are, split
        ! them up into parts such that each face has its own unique vertex pair.

        ! Notes
        !======
        ! Note 1: it is assumed that saddle points which should be treated with
        ! equal field value have exactly the same field value.

        ! Note 2: for tangency points, we need to be careful when tracing contours.
        ! The initial or end part of the contour, which should start or end at the
        ! tangency point, should lie strictly inside the vessel, i.e. the first
        ! point should be inside the vessel domain. If this is not the case, these
        ! points are eliminated up to the first point that lies inside of the
        ! vessel again. Problems may arise if contour mesh resolution is (locally)
        ! not fine enough, and mesh refinement may be needed. However, this may
        ! require a very fine mesh in cases where the contour and vessel are almost
        ! aligned over a substantial part of the domain. Therefore, we propose a
        ! different approach: first, we compute the contours. Then, we check which
        ! boundary faces (type 3) the contours intersect. We then process the
        ! contours and boundary faces as follows: only segments are kept that start
        ! from the tangency point. However, if the first of such segments goes
        ! outside of the vessel, we have to treat it differently (this can happen
        ! due to different order of approximation of the contour (first order) and
        ! determination of tangency points (higher order)). Then, we remove this
        ! first segment from tangency point to the intersection and connect the
        ! tangency point to the first next point of the contour. We do the same for
        ! the vessel polygon, which ensures a proper and expected topology.

        ! Note 3: other boundaries and points may be inserted afterwards that still
        ! lie outside of the domain (this is typically the case with saddle point
        ! contours). These contours are not trimmed beforehand since we still need
        ! to find the proper intersections with the vessel wall. 

        ! Note 4: also for saddle points/cut-like boundaries, we need
        ! to be careful when computing intersections. It is possible that
        ! near the x-point, the separatrix intersects with a core boundary
        ! due to inaccurate tracing of the field line contour locally. 
        ! This needs to be hedged for by modifying the contours/lines 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        real(R8)                                :: startsr, endsr
        real(R8), allocatable, dimension(:)     :: pspx, pspy, pspf, &
            xout, yout, iout, jout, tscr, tsfr, tx, ty, tsrfh, tsrsh
        integer(I8)                             :: npsp, &
            intersectind, nstc, ntpc, nint, startind, endind, indtpc
        integer(I8), allocatable, dimension(:)  :: psptype, &
            sortind, pspID, allcurvetypes, allfsIDs, &
            vindI, vindJ, notdelind, vindIfh, vindIsh, tsc, tfaceind, &
            tsf, tcontourind, tcstartind, tcendind
        logical                                 :: intfacestart, &
            intcstart, dointersect
        logical, allocatable, dimension(:)      :: tracepoints, &
            ispseudosaddlepoint, isnbface, rmind, isstartingcontour, &
            isendingcontour, forcevert

        type(ContourUDT), allocatable           :: tc(:), allc(:), alltpc(:)
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        type(IntegerDynamicArrayUDT)            :: curvetypes, fsIDs
        type(IntegerDynamicArrayUDT), allocatable       :: sc(:), &
            faceind(:), sf(:), contourind(:)
        type(RealDynamicArrayUDT), allocatable          :: scr(:), sfr(:)

        type(PolygonSetUDT)                             :: bndps
        type(PolygonUDT), allocatable, dimension(:)     :: bndpol
        class(PLF2DOptionsUDT), allocatable     :: bndplfoptions

        ! Loop 
        integer(I8)                             :: i, j, k, cc
        
        ! Initialize
        !===========
        ! Associate
        associate(nfs       => topomesh%nfs)
            
        ! Initialize the contour structure & dynamic arrays
        allocate(allc(0))
        curvetypes = ConstructIntegerDynamicArray()
        fsIDs = ConstructIntegerDynamicArray()

        ! Build 'saddle' points
        ispseudosaddlepoint = (topomesh%vert%type == TMvertexmaxID) .or. &
            (topomesh%vert%type == TMvertexminID) .or. &
            (topomesh%vert%type == TMvertexsaddleID) .or. & 
            (topomesh%vert%type == TMvertextp1ID) .or. & 
            (topomesh%vert%type == TMvertextp2ID) .or. &
            (topomesh%vert%type == TMvertexsplitID)
        npsp = count(ispseudosaddlepoint)
        allocate(pspx(npsp), pspy(npsp), pspf(npsp), psptype(npsp), &
            pspID(npsp))
        pspx    = pack(topomesh%vert%x, ispseudosaddlepoint)
        pspy    = pack(topomesh%vert%y, ispseudosaddlepoint) 
        pspf    = pack(topomesh%vert%fval, ispseudosaddlepoint) 
        psptype = pack(topomesh%vert%type, ispseudosaddlepoint)
        pspID   = pack(topomesh%vert%ID, ispseudosaddlepoint)

        ! Check for saddle points that are connected - don't trace twice
        allocate(tracepoints(npsp))
        tracepoints = .true.
        !do i = 1, npsp-1
        !    do j = i+1, npsp
        !        if ((topomesh%vert%fsID(i) == topomesh%vert%fsID(j)) .and. &
        !            (psptype(i) == TMvertexsaddleID) .and. (psptype(j) == TMvertexsaddleID)) then 
        !            tracepoints(j) = .false.
        !        end if 
        !    end do 
        !end do

        ! Add saddle points to field tracer - only tangency points and
        ! extrema, no other 'regular' vertices

        fieldtracer%xs = pspx 
        fieldtracer%ys = pspy
        fieldtracer%vs = pspf
        fieldtracer%IDs = pspID
        fieldtracer%order = 0*psptype  
        fieldtracer%issaddle = topomesh%vert%type(pspID) == TMvertexsaddleID 

        ! Trace contours
        !===============
        ! Add saddle point contours
        !-------------------------- 
        do i = 1, npsp
            if ((psptype(i) == TMvertexsaddleID) .and. tracepoints(i)) then 
                ! Trace
                tc = fieldtracer%TraceContours([pspx(i)], [pspy(i)])
                
                ! Process
                call CleanContours(tc)

                ! Check if any other saddle points were encountered 
                ! during tracing. If so, do not trace these anymore 
                ! (leads to duplicate faces)
                do j = 1, size(tc)
                    if (tc(j)%startsaddle /= 0) then
                        if (psptype(tc(j)%startsaddle) == TMvertexsaddleID) then  
                            tracepoints(tc(j)%startsaddle) = .false. 
                        end if 
                    end if 
                    if (tc(j)%endsaddle /= 0) then
                        if (psptype(tc(j)%endsaddle) == TMvertexsaddleID) then  
                            tracepoints(tc(j)%endsaddle) = .false. 
                        end if 
                    end if 
                end do 

                ! Check intersections with radial core boundaries
                do j = 1, topomesh%face%ntot
                    ! Check where exactly it (should) intersect
                    if ((topomesh%face%type(j) == TMfaceradID)) then 
                        
                        ! Loop over all contours
                        do k = 1, size(tc)
                            ! Check for intersections
                            dointersect = .true.
                            if ((topomesh%face%vert(j, 1) == tc(k)%startsaddle)) then
                                intfacestart = .true.
                                intcstart = .true.
                            elseif ((topomesh%face%vert(j, 2) == tc(k)%startsaddle)) then
                                intfacestart = .false.
                                intcstart = .true.
                            elseif ((topomesh%face%vert(j, 1) == tc(k)%endsaddle)) then
                                intfacestart = .true.
                                intcstart = .false.
                            elseif ((topomesh%face%vert(j, 2) == tc(k)%endsaddle)) then
                                intfacestart = .false.
                                intcstart = .false.
                            else 
                                ! No intersection found - cycle
                                dointersect = .false. 
                            end if 

                            if (dointersect) then 
                                call SimplePolygonIntersections(tc(k)%x, tc(k)%y, &
                                    topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), &
                                    xout, yout, vindI, vindJ, iout, jout)

                                ! Check intersections - first one should be in x-point itself
                                ! Note: we don't need to hedge for closed contours here yet
                                nstc = size(tc(k)%x)-1 ! number of segments
                                if (allocated(iout)) then 
                                    if (size(iout) < 1) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'could not find any intersection ' // & 
                                            'between radial line and contour ' // & 
                                            'though they should start or end ' // &
                                            'in the same point')
                                    end if 
                                else 
                                    call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'could not find any intersection ' // & 
                                            'between radial line and contour ' // & 
                                            'though they should start or end ' // &
                                            'in the same point')
                                end if 

                                ! If any intersections are found, check
                                if (intfacestart .and. intcstart) then 
                                    if (.not. (jout(1) == 0.0_R8) .and. .not. (iout(1) == 0.0_R8)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (intfacestart .and. .not. intcstart) then 
                                    if (.not. (jout(1) == 0.0_R8) .and. .not. (iout(size(iout)) == nstc)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (.not. intfacestart .and. intcstart) then 
                                    if (.not. (jout(size(jout)) == topomesh%face%x(j)%Size()-1) .and. &
                                        .not. (iout(1) == 0.0_R8)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if 
                                if (.not. intfacestart .and. .not. intcstart) then 
                                    if (.not. (jout(size(jout)) == topomesh%face%x(j)%Size()-1) .and. &
                                        .not. (iout(size(iout)) == nstc)) then 
                                        call gdErrorHandler('AddTopologicalMeshContour: ' // & 
                                            'first intersection between separatrix contour ' // & 
                                            'and radial line should be in first ' // & 
                                            'point of separatrix but this is not the case')
                                    end if
                                end if

                                ! Sort the intersections according to contour
                                sortind = [(k, k = 1, size(xout))]
                                call Sort(iout, ind=sortind, ascend=.true.)
                                vindI = vindI(sortind)
                                
                                ! Sort the intersections according to face
                                sortind = [(k, k = 1, size(xout))]
                                call Sort(jout, ind=sortind, ascend=.true.)
                                vindJ = vindJ(sortind)

                                ! If there are other intersections, adjust the face and the contour
                                if (size(xout) > 1) then 
                                    ! Face
                                    !-----
                                    tx = topomesh%face%x(j)%Get()
                                    ty = topomesh%face%y(j)%Get()
                                    if (intfacestart) then 
                                        ! Face is oriented from start to end
                                        call DeleteCurveSegment(tx, ty, [maxval(jout)], 'start', &
                                            [distfrac, distfrac], .true., .true.)
                                        ! delind = [(cc, cc = 2, maxval(vindJ))]
                                    else
                                        call DeleteCurveSegment(tx, ty, [jout(1)], 'end', [distfrac, distfrac], .true., .true.)
                                        ! delind = [(cc, cc = vindJ(1)+1, topomesh%face%x(j)%Size()-1)]
                                    end if 

                                    ! Delete
                                    !call topomesh%face%x(j)%Remove(delind)
                                    !call topomesh%face%y(j)%Remove(delind)
                                    call topomesh%face%x(j)%Set(tx)
                                    call topomesh%face%y(j)%Set(ty)

                                    ! Reconvert to polygon
                                    allocate(forcevert(size(tx)))
                                    forcevert([1, size(tx)]) = .true. ! force final vertices to remain the same
                                    call topomesh%face%pol(j)%Construct(&
                                        topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), forcevert_opt=forcevert)
                                    deallocate(forcevert)

                                    ! Contour
                                    !--------
                                    ! Need to hedge for closed contour
                                    if (tc(k)%isclosed) then
                                        ! Sanity check
                                        if ((iout(1) /= 0.0_R8) .or. &
                                            (iout(size(iout)) /= size(tc(k)%x)-1)) then 
                                            ! Intersection should happen in both end and start point
                                            print *, 'vertex : ', i
                                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                                'closed contour does not intersect in both ' // &
                                                'start and en points, check input')
                                        end if

                                        ! Check if there are exactly two
                                        ! intersections - in that case, 
                                        ! this is simply a radial line
                                        ! intersecting with a closed 
                                        ! contour
                                        if (size(iout) /= 2) then 
                                            ! If not, issue warning
                                            print *, 'vertex: ', i 
                                            print *, 'AddTopologicalMeshContours: ' // & 
                                                'multiple intersections with ' // & 
                                                'closed separatrix found, attempting ' // & 
                                                'to adjust separatrix and radial line' 
                                        end if 

                                        

                                        ! Take only intersections that 
                                        ! are not in end points
                                        vindI = vindI(2:size(vindI)-1)

                                        ! Partition into intersections in 
                                        ! first half of curve and second half
                                        allocate(vindIfh(count(vindI < nstc/2_I8)), vindIsh(count(vindI >= nstc/2_I8)))
                                        vindIfh = pack(vindI, vindI < nstc/2_I8)
                                        vindIsh = pack(vindI, vindI >= nstc/2_I8)

                                        ! Check which part to keep
                                        if ((size(vindIfh) > 0) .and. (size(vindIsh) > 0)) then 
                                            ! Intersections at both sides, take 'middle' piece
                                            notdelind = [1, (cc, cc = maxval(vindIfh)+1, minval(vindIsh)-1), nstc+1]
                                        elseif (size(vindIfh) > 0) then 
                                            ! Intersection only at first half
                                            notdelind = [1, (cc, cc = maxval(vindIfh)+1, nstc+1)]
                                        elseif (size(vindIsh) > 0) then 
                                            ! Intersection only at second half
                                            notdelind = [(cc, cc = 1, minval(vindIsh)-1), nstc+1]
                                        else
                                            ! Only two intersections in 
                                            ! start and end, move along
                                            notdelind = [(cc, cc = 1, size(tc(k)%x))]

                                        end if 
                                        deallocate(vindIfh, vindIsh)

                                    elseif (intcstart) then 
                                        ! Contour is oriented from start to end
                                        notdelind = [1, (cc, cc = maxval(vindI)+1, size(tc(k)%x))]
                                    else
                                        notdelind = [(cc, cc = 1, vindI(1)), size(tc(k)%x)]
                                    end if  

                                    ! Remove
                                    tc(k)%x = tc(k)%x(notdelind)
                                    tc(k)%y = tc(k)%y(notdelind)
                                end if 
                            end if 
                        end do 
                    end if 
                end do

                ! Add
                allc = [allc, tc]
                call curvetypes%Append(spread(TMfacesepID, 1, size(tc)))
                
                ! Add flux surface ID - should be X-point flux surface ID 
                ! nfs = nfs + 1
                call fsIDs%Append(spread(topomesh%vert%fsID(pspID(i)), 1, size(tc)))
            end if 
        end do

        ! Add tangency point contours
        !----------------------------
        ! First, trace all contours
        allocate(alltpc(0))
        do i = 1, npsp
            if ((psptype(i) == TMvertextp2ID) .and. tracepoints(i)) then 
                ! Trace contour
                tc = fieldtracer%TraceContours([pspx(i)], [pspy(i)])

                ! Process
                call CleanContours(tc)

                ! Check if the tangency point is connected to any other
                ! points (we include the current point because that doesn't matter)
                ! These points should not be traced anymore because the
                ! contour algorithm should've gone through all of them
                do k = 1, size(tc)
                    if (tc(k)%startsaddle /= 0) then 
                        tracepoints(tc(k)%startsaddle) = .false. 
                    end if 
                    if (tc(k)%endsaddle /= 0) then 
                        tracepoints(tc(k)%endsaddle) = .false. 
                    end if 
                end do 

                ! Sanity checks
                do k = 1, size(tc)
                    ! Is at least one of both end points a saddle point? 
                    ! If not, throw an error, this shouldn't happen
                    if (tc(k)%startsaddle == 0 .and. tc(k)%endsaddle == 0) then 
                        call WriteTopologicalMesh(topomesh, 'topomesh_error')
                        call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                            'tangency point contour detected that does not start ' // & 
                            'or end in any topomesh vertex, this is likely a bug. ' // & 
                            'Current topomesh written in topomesh_error.dat')
                    end if 
                end do

                ! Add
                alltpc = [alltpc, tc]
                call curvetypes%Append(spread(TMfacepolID, 1, size(tc)))
                
                ! Add flux surface ID
                call fsIDs%Append(spread(i, 1, size(tc)))
            end if 
        end do 

        ! Compute all intersections with other faces
        ntpc = size(alltpc)
        allocate(sc(ntpc), scr(ntpc), faceind(ntpc), sf(topomesh%face%ntot), &
            sfr(topomesh%face%ntot), contourind(topomesh%face%ntot))
        do i = 1, ntpc
            sc(i) = ConstructIntegerDynamicArray()
            scr(i) = ConstructRealDynamicArray()
            faceind(i) = ConstructIntegerDynamicArray()
        end do 
        do i = 1, topomesh%face%ntot
            sf(i) = ConstructIntegerDynamicArray()
            sfr(i) = ConstructRealDynamicArray()
            contourind(i) = ConstructIntegerDynamicArray()
        end do 

        ! Loop
        !$omp parallel do default(none) &
        !$omp schedule(dynamic) collapse(2) & 
        !$omp private(xout, yout, vindI, vindJ, iout, jout) & 
        !$omp shared(ntpc, topomesh, alltpc, sc, scr, faceind, sf, sfr, contourind)
        do i = 1, ntpc
            do j = 1, topomesh%face%ntot
                if (topomesh%face%type(j) == TMfacebndID) then 
                    ! Compute intersections
                    call SimplePolygonIntersections(alltpc(i)%x, alltpc(i)%y, &
                        topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), &
                        xout, yout, vindI, vindJ, iout, jout)

                    ! Add, if any
                    if (allocated(xout)) then 
                        if (size(xout) > 0) then 
                            !$omp critical
                            ! Add intersection data to contour
                            call sc(i)%Append(vindI)
                            call scr(i)%Append(iout)
                            call faceind(i)%Append(spread(j, 1, size(vindI)))

                            ! Add intersection data to face
                            call sf(j)%Append(vindJ)
                            call sfr(j)%Append(jout)
                            call contourind(j)%Append(spread(i, 1, size(vindI))) 
                            !$omp end critical
                        end if 
                    end if
                end if
            end do 
        end do
        !$omp end parallel do 

        ! Process contours
        do i = 1, ntpc
            ! Unpack
            tsc = sc(i)%Get()
            tfaceind = faceind(i)%Get()
            tscr = scr(i)%Get()
            nint = size(tsc)
            nstc = size(alltpc(i)%x)-1

            ! Sort 
            sortind = [(k, k = 1, nint)]
            call Sort(tscr, ind=sortind, ascend=.true.)
            tsc = tsc(sortind)
            tfaceind = tfaceind(sortind)

            ! Check
            if (tscr(1) /= 0.0_R8) then 
                ! First intersection should always be in tangency point
                call gdErrorHandler('AddTopologicalMeshContours: first ' // & 
                    'intersection of tangency point contour is not in  ' // & 
                    'tangency point itself' )
            end if 
            if (alltpc(i)%isclosed) then 
                ! Last intersection should also be in tangency point
                if (tscr(nint) /= size(alltpc(i)%x)-1) then 
                    call gdErrorHandler('AddTopologicalMeshContours: last ' // & 
                    'intersection of closed tangency point contour is not in  ' // & 
                    'tangency point itself' )
                end if 
            end if 

            ! Determine which intersections were in a 'neighbouring' face
            ! (i.e. a boundary face that has the tangency point)
            allocate(isnbface(nint))
            isnbface = .false.
            do j = 1, nint
                if (any(topomesh%face%vert(tfaceind(j), :) == alltpc(i)%startsaddle)) then 
                    isnbface(j) = .true.
                end if 
                if (any(topomesh%face%vert(tfaceind(j), :) == alltpc(i)%endsaddle) .and. &
                    alltpc(i)%endsaddle /= 0) then 
                    isnbface(j) = .true.
                end if 
            end do             

            ! Determine which part(s) of contour(s) to keep
            if (alltpc(i)%isclosed) then ! closed contour
                ! This is much more tricky, need to check additional cases
                if (all(isnbface)) then 
                    ! No intersections with other boundaries
                    if (nint == 2) then 
                        ! Simply in start and end - add full contour, so 
                        ! do nothing
                    else
                        ! Here, we need to look at which face IDs the 
                        ! intersections have, since it is 
                        ! Here, we don't have much more to go on than 
                        ! assuming that these intersections are somewhere
                        ! in the start/end nodes of the contour and that
                        ! we can partition based on that. This is of course
                        ! not a general way and may fail if the amount of
                        ! contour points is too low
                        ! Print a warning
                        print *, 'AddTopologicalMeshContours: multiple ' // & 
                            'intersections found for closed contour with ' // & 
                            'only neighbouring boundary faces, attempting to ' // &
                            'split contour by partitioning intersection into ' // & 
                            'intersections at start and end. This may not result ' // & 
                            'in desired behavior... (vertex: ', alltpc(i)%startsaddle, ')'

                        ! Determine start & end intersections
                        allocate(vindIfh(count(tsc < nstc/2)), &
                            vindIsh(count(tsc >= nstc/2)))
                        allocate(tsrfh(size(vindIfh)), tsrsh(size(vindIsh)))
                        vindIfh = pack(tsc, tsc < nstc/2)
                        vindIsh = pack(tsc, tsc >= nstc/2)
                        tsrfh = pack(tscr, tsc < nstc/2)
                        tsrsh = pack(tscr, tsc >= nstc/2)

                        ! Determine start index
                        if (size(vindIfh) > 0) then 
                            startind = maxval(vindIfh)
                            startsr = maxval(tsrfh)
                        else
                            startind = 2
                            startsr = 0.0_R8
                        end if 

                        ! Determine end index
                        if (size(vindIsh) > 0) then 
                            endind = minval(vindIsh)+1
                            endsr = minval(tsrsh)
                        else
                            endind = nstc
                            endsr = real(nstc, kind=R8)
                        end if 
                        
                        ! Adjust contour
                        if (startsr /= 0.0_R8 .and. endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr, endsr], 'both', [distfrac, distfrac], .true., .true.)
                        elseif (startsr /= 0.0_R8) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr], 'start', [distfrac, distfrac], .true., .true.)
                        elseif (endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [endsr], 'end', [distfrac, distfrac], .true., .true.)
                        end if 
                        !alltpc(i)%x = alltpc(i)%x([1, (k, k = startind, endind), nstc+1])
                        !alltpc(i)%y = alltpc(i)%y([1, (k, k = startind, endind), nstc+1])

                        ! Housekeeping
                        deallocate(vindIfh, vindIsh, tsrfh, tsrsh)
                    end if
                else
                    ! At least one intersection with another boundary. 
                    ! Note: here we do want to keep the intersection 
                    ! with these other boundaries in the contour!
                    ! Need to find first and second segment 

                    ! First segment
                    endind = findloc(isnbface, .false., 1, back=.false.)

                    ! Add this segment as additional contour
                    alltpc = [alltpc, alltpc(i)]
                    indtpc = size(alltpc)
                    alltpc(indtpc)%isclosed = .false.
                    alltpc(indtpc)%endsaddle = 0 ! doesn't end anymore in saddle point
                    if (tscr(endind-1) == 0.0_R8) then 
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind)], 'end', [-distfrac], .true., .false.)
                    else
                        ! Also need to delete a first part
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind-1:endind)], 'both', [distfrac, -distfrac], .true., .false.)
                    end if 
                    !call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                    !    [tscr(endind)], 'end', [0.0_R8], .true., .false.)
                    !alltpc(indtpc)%x = alltpc(i)%x([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])
                    !alltpc(indtpc)%y = alltpc(i)%y([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])

                    ! Append flux surface ID etc as well!
                    call curvetypes%Append(curvetypes%Get(size(allc) + i))
                    call fsIDs%Append(fsIDs%Get(size(allc) + i))

                    ! Second segment
                    startind = findloc(isnbface, .false., 1, back=.true.)

                    ! Add this segment by adjusting existing contour
                    alltpc(i)%isclosed = .false.
                    alltpc(i)%startsaddle = 0 ! doesn't start anymore in saddle point
                    if (tscr(startind+1) == real(nstc, kind=R8)) then 
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind)], 'start', [-distfrac], .false., .true.)
                    else
                        ! Also need to delete last part
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind:startind+1)], 'both', [-distfrac, distfrac], .false., .true.)
                    end if
                    !alltpc(i)%x = alltpc(i)%x([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    !alltpc(i)%y = alltpc(i)%y([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    
                end if 
            elseif (alltpc(i)%startsaddle /= 0 .and. alltpc(i)%endsaddle /= 0) then  ! open contour segment but with two tangency points
                if (all(isnbface)) then 
                    ! No intersections with other boundaries
                    if (nint == 2) then 
                        ! Simply in start and end - add full contour, so 
                        ! do nothing
                    else
                        ! Here, we don't have much more to go on than 
                        ! assuming that these intersections are somewhere
                        ! in the start/end nodes of the contour and that
                        ! we can partition based on that. This is of course
                        ! not a general way and may fail if the amount of
                        ! contour points is too low
                        ! Print a warning
                        print *, 'AddTopologicalMeshContours: multiple ' // & 
                            'intersections found for closed contour with ' // & 
                            'only neighbouring boundary faces, attempting to ' // &
                            'split contour by partitioning intersection into ' // & 
                            'intersections at start and end. This may not result ' // & 
                            'in desired behavior... (vertex: ', alltpc(i)%startsaddle, ')'

                        ! Determine start & end intersections
                        allocate(vindIfh(count(tsc < nstc/2)), &
                            vindIsh(count(tsc >= nstc/2)))
                        allocate(tsrfh(size(vindIfh)), tsrsh(size(vindIsh)))
                        vindIfh = pack(tsc, tsc < nstc/2)
                        vindIsh = pack(tsc, tsc >= nstc/2)
                        tsrfh = pack(tscr, tsc < nstc/2)
                        tsrsh = pack(tscr, tsc >= nstc/2)

                        ! Determine start index
                        if (size(vindIfh) > 0) then 
                            startind = maxval(vindIfh)
                            startsr = maxval(tsrfh)
                        else
                            startind = 2
                            startsr = 0.0_R8
                        end if 

                        ! Determine end index
                        if (size(vindIsh) > 0) then 
                            endind = minval(vindIsh)+1
                            endsr = minval(tsrsh)
                        else
                            endind = nstc
                            endsr = real(nstc, kind=R8)
                        end if 
                        
                        ! Adjust contour
                        if (startsr /= 0.0_R8 .and. endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr, endsr], 'both', [distfrac, distfrac], .true., .true.)
                        elseif (startsr /= 0.0_R8) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [startsr], 'start', [distfrac, distfrac], .true., .true.)
                        elseif (endsr /= real(nstc, kind=R8)) then 
                            call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                                [endsr], 'end', [distfrac, distfrac], .true., .true.)
                        end if 
                        !alltpc(i)%x = alltpc(i)%x([1, (k, k = startind, endind), nstc+1])
                        !alltpc(i)%y = alltpc(i)%y([1, (k, k = startind, endind), nstc+1])

                        ! Housekeeping
                        deallocate(vindIfh, vindIsh, tsrfh, tsrsh)
                    end if
                else
                    ! At least one intersection with another boundary. 
                    ! Note: here we do want to keep the intersection 
                    ! with these other boundaries in the contour!
                    ! Need to find first and second segment 

                    ! First segment
                    endind = findloc(isnbface, .false., 1, back=.false.)

                    ! Add this segment as additional contour
                    alltpc = [alltpc, alltpc(i)]
                    indtpc = size(alltpc)
                    alltpc(indtpc)%isclosed = .false.
                    alltpc(indtpc)%endsaddle = 0 ! doesn't end anymore in saddle point
                    if (tscr(endind-1) == 0.0_R8) then 
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind)], 'end', [-distfrac], .true., .false.)
                    else
                        ! Also need to delete a first part
                        call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                            [tscr(endind-1:endind)], 'both', [distfrac, -distfrac], .true., .false.)
                    end if 
                    !call DeleteCurveSegment(alltpc(indtpc)%x, alltpc(indtpc)%y, &
                    !    [tscr(endind)], 'end', [0.0_R8], .true., .false.)
                    !alltpc(indtpc)%x = alltpc(i)%x([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])
                    !alltpc(indtpc)%y = alltpc(i)%y([1, (k, k = tsc(startind-1)+1, tsc(startind)+1)])

                    ! Append flux surface ID etc as well!
                    call curvetypes%Append(curvetypes%Get(size(allc) + i))
                    call fsIDs%Append(fsIDs%Get(size(allc) + i))

                    ! Second segment
                    startind = findloc(isnbface, .false., 1, back=.true.)

                    ! Add this segment by adjusting existing contour
                    alltpc(i)%isclosed = .false.
                    alltpc(i)%startsaddle = 0 ! doesn't start anymore in saddle point
                    if (tscr(startind+1) == real(nstc, kind=R8)) then 
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind)], 'start', [-distfrac], .false., .true.)
                    else
                        ! Also need to delete last part
                        call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                            [tscr(startind:startind+1)], 'both', [-distfrac, distfrac], .false., .true.)
                    end if
                    !alltpc(i)%x = alltpc(i)%x([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    !alltpc(i)%y = alltpc(i)%y([(k, k = tsc(endind), tsc(endind+1)), nstc+1])
                    
                end if 

            else ! open contour with only one tangency point that is located at the start
                ! Check which intersection is the last intersection with
                ! the neighbouring boundary (should be first one)
                intersectind = findloc(isnbface, .true., 1, back=.true.)

                ! Sanity checks
                if (intersectind == 0) then 
                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        'tangency point contour intersects in tangency point ' // & 
                        'but not in neighbouring face - this is likely a bug')
                end if 
                if (any(.not. isnbface(1:intersectind))) then 
                    ! Issue warning - we got intersections of the contour
                    ! with non-neighbour faces inbetween - this is unexpected
                    print *, 'AddTopologicalMeshContours: intersections ' // & 
                        'found for tangency point ', alltpc(i)%startsaddle, &
                        'that occur inbetween intersections with contour and ' // & 
                        'neighbouring boundary face. Unexpected, intersections ' // & 
                        'are removed'
                end if 
                if (intersectind == size(isnbface)) then 
                    ! Unexpected, only intersections in neighbouring faces 
                    ! but not in other boundary face
                    call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        'tangency point contour intersects in tangency point ' // & 
                        'and neighbouring face but not in other boundary ' // &
                        'faces - unexpected')
                end if 

                ! Keep only this part of the contour coordinates
                call DeleteCurveSegment(alltpc(i)%x, alltpc(i)%y, &
                    [tscr(intersectind:intersectind+1)], 'both', [distfrac, -distfrac], .true., .false.)
                !notdelind = [(k, k = 2, tsc(intersectind))]
                !alltpc(i)%x = alltpc(i)%x(notdelind)
                !alltpc(i)%y = alltpc(i)%y(notdelind)
            end if 

            ! Housekeeping
            deallocate(isnbface)

        end do 

        ! Process faces
        do i = 1, topomesh%face%ntot 
            if (topomesh%face%type(i) == TMfacebndID) then 
                ! Unpack
                tsf = sf(i)%Get()
                tsfr = sfr(i)%Get()
                tcontourind = contourind(i)%Get()
                nint = size(tsf)

                ! Sort
                sortind = [(k, k = 1, nint)]
                call Sort(tsfr, ind=sortind, ascend=.true.)
                tsf = tsf(sortind)
                tcontourind = tcontourind(sortind)

                ! Check intersections
                if (nint > 0) then 
                    ! Check which intersections should be removed
                    allocate(rmind(nint))
                    rmind = .false.

                    do j = 1, nint
                        if (&
                            any(topomesh%face%vert(i, :) == alltpc(tcontourind(j))%startsaddle) .or. &
                            any(topomesh%face%vert(i, :) == alltpc(tcontourind(j))%endsaddle)) then 
                            rmind(j) = .true.
                        end if
                    end do 

                    ! Check how to remove intersections 
                    if (all(topomesh%vert%type(topomesh%face%vert(i, : )) == TMvertextp2ID)) then 
                        !print *, 'face ID: ', i
                        !print *, 'face vertices: ', topomesh%face%vert(i, :) 
                        !call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                        !    'boundary detected with two type 2 tangency points, ' // & 
                        !    'this is likely a bug')

                        ! Sanity checks
                        if ((tsfr(1) /= 0.0_R8) .or. (tsfr(nint) /= topomesh%face%x(i)%Size()-1)) then 
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'at start and end, but not all intersections found')
                        end if 

                        ! Check which contours have the start and end
                        ! vertex as a starting/ending point
                        isstartingcontour = tsfr == 0.0_R8
                        isendingcontour = tsfr == topomesh%face%x(i)%Size()-1
                        allocate(tcstartind(count(isstartingcontour)), tcendind(count(isendingcontour)))
                        tcstartind = pack(tcontourind, isstartingcontour)
                        tcendind = pack(tcontourind, isendingcontour)

                        ! Check which intersections these are, these
                        ! give limits for face trimming
                        startind = 0
                        endind = nint+1 
                        k = 1
                        do while (k <= nint)
                            if (.not. any(tcontourind(k) == tcstartind)) then 
                                ! Exit the loop
                                exit 
                            else
                                ! Update counters
                                k = k + 1
                                startind = startind + 1
                            end if 
                        end do 
                        k = nint 
                        do while (k > 0)
                            if (.not. any(tcontourind(k) == tcendind)) then 
                                ! Exit the loop
                                exit
                            else
                                ! Update counters
                                k = k - 1
                                endind = endind - 1
                            end if 
                        end do 

                        ! Trim face
                        tx = topomesh%face%x(i)%Get()
                        ty = topomesh%face%y(i)%Get()
                        if (tsfr(startind) /= 0.0_R8 .and. tsfr(endind) /= real(topomesh%face%x(i)%Size()-1)) then 
                            call DeleteCurveSegment(tx, ty, tsfr([startind, endind]), &
                                'both', [distfrac, distfrac], .true., .true.)
                        elseif (tsfr(startind) /= 0.0_R8) then 
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'start', [distfrac, distfrac], .true., .true.)
                        elseif (tsfr(endind) /= real(topomesh%face%x(i)%Size()-1, kind=R8)) then 
                            call DeleteCurveSegment(tx, ty, tsfr([endind]), &
                                'end', [distfrac, distfrac], .true., .true.)
                        end if 

                            ! Reconstruct polygon
                        call topomesh%face%pol(i)%Construct(&
                            topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        call topomesh%face%x(i)%Set(tx)
                        call topomesh%face%y(i)%Set(ty)
                        !call topomesh%face%x(i)%Remove([(k, k = 2, tsf(startind)+1), &
                        !    (k, k = tsf(endind+1), topomesh%face%x(i)%Size()-1)])
                        !call topomesh%face%y(i)%Remove([(k, k = 2, tsf(startind)+1), &
                        !    (k, k = tsf(endind+1), topomesh%face%y(i)%Size()-1)])

                        ! Reconstruct polygon
                        call topomesh%face%pol(i)%Construct(&
                            topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())

                        ! Housekeeping
                        deallocate(tcstartind, tcendind)

                    elseif (topomesh%vert%type(topomesh%face%vert(i, 1)) == TMvertextp2ID) then

                        ! Sanity checks
                        if (tsfr(1) /= 0.0_R8) then
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'initially, but no intersection found')
                        end if 

                        ! Find location
                        startind = findloc(rmind, .true., 1, back=.true.)

                        ! Sanity check
                        if (startind == 0) then 
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                'face has intersections but none to be removed - unexpected')
                        end if 

                        ! Trim face
                        if (tsfr(startind) /= 0.0_R8) then 
                            tx = topomesh%face%x(i)%Get()
                            ty = topomesh%face%y(i)%Get()
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'start', [distfrac, distfrac], .true., .true.)
                            call topomesh%face%x(i)%Set(tx)
                            call topomesh%face%y(i)%Set(ty)

                            ! Reconstruct polygon
                            call topomesh%face%pol(i)%Construct(&
                                topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        end if 

                        !call topomesh%face%x(i)%Remove([(k, k = 2, tsf(startind)+1)])
                        !call topomesh%face%y(i)%Remove([(k, k = 2, tsf(startind)+1)])

                        

                    elseif (topomesh%vert%type(topomesh%face%vert(i, 2)) == TMvertextp2ID) then

                        ! Sanity checks
                        if (tsfr(nint) /= topomesh%face%x(i)%Size()-1) then
                            print *, 'face ID: ', i
                            print *, 'face vertices: ', topomesh%face%vert(i, :) 
                            call gdErrorHandler('AddTopologicalMeshContours: ' //& 
                                'face should intersect with tangency point ' // & 
                                'initially, but no intersection found')
                        end if 

                        ! Find location
                        startind = findloc(rmind, .true., 1, back=.false.)

                        ! Sanity check
                        if (startind == 0) then 
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshContours: ' // & 
                                'face has intersections but none to be removed - unexpected')
                        end if 

                        ! Trim face
                        if (tsfr(startind) < real(topomesh%face%x(i)%Size()-1, kind=R8)) then 
                            tx = topomesh%face%x(i)%Get()
                            ty = topomesh%face%y(i)%Get()
                            call DeleteCurveSegment(tx, ty, tsfr([startind]), &
                                'end', [distfrac, distfrac], .true., .true.)
                            call topomesh%face%x(i)%Set(tx)
                            call topomesh%face%y(i)%Set(ty)

                            ! Reconstruct polygon
                            call topomesh%face%pol(i)%Construct(&
                                topomesh%face%x(i)%Get(), topomesh%face%y(i)%Get())
                        end if 
                        !call topomesh%face%x(i)%Remove([(k, k = tsf(startind), topomesh%face%x(i)%Size()-1)])
                        !call topomesh%face%y(i)%Remove([(k, k = tsf(startind), topomesh%face%y(i)%Size()-1)])

                        
                    else
                        ! This can happen, just skip since we shouldn't 
                        ! adjust anything here
                    end if

                    ! Housekeeping
                    deallocate(rmind)

                end if 

                
            end if 
        end do 

        ! Add contours (IDs etc should be done already...)
        allc = [allc, alltpc]

        


        ! Add to the topology mesh
        !=========================
        ! Clean the contours (just to be sure)
        call CleanContours(allc)

        if (options%writedebugoutput) then 
            allocate(allpc(size(allc)))
            do i = 1, size(allc)
                call allpc(i)%Construct(allc(i)%x, allc(i)%y)
            end do 
            call tempps%Construct(allpc)
            call tempps%WriteData('topocontours_all_beforeinsertion')
        end if

        ! Add the contours
        allcurvetypes = curvetypes%Get()
        allfsIDs = fsIDs%Get()
        do i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), allcurvetypes(i), allfsIDs(i))
        end do 

        ! Update the vessel description again
        allocate(bndpol(count(topomesh%face%type == TMfacebndID)))
        bndpol = pack(topomesh%face%pol, topomesh%face%type == TMfacebndID)
        call bndps%Construct(bndpol)
        call ConstructVesselPolygonSet(vessel, bndps, .true.)
        allocate(PLF2DClosedExactOptionsUDT::bndplfoptions)
        call InitializePolygonLevelsetFunction2D(vessel%plfvessel, vessel%polygonset, bndplfoptions)
        deallocate(bndpol)

        ! Trim the topological mesh
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_temp')
        end if 
        call TrimTopologicalMesh(topomesh, magneticField, vessel)
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Split boundaries
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_temp')
        end if 
        call SplitTopologicalMeshFaces(topomesh)   
        
        ! Boundary split vertex contours
        !===============================
        call AddBoundarySplitVertexContours(topomesh, magneticField, vessel, &
            fieldtracer)

        call SimplifyTopologicalMeshFaces(topomesh)
        
        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Contour insertion into topological mesh
    subroutine InsertTopologicalMeshContour(topomesh, magneticField, contour, &
        contourtype, contourfsID, remfout) 

        ! Description
        !============
        ! This routine provides a general way to insert a curve into an existing
        ! topological mesh. The topological mesh should contain vertices and faces
        ! already, but no cells yet (these are also not updated here). The curve
        ! should have at least x and y data and, if available, start and/or end
        ! vertex IDs of vertices that occur in the topological mesh. If new
        ! vertices should be introduced, this should be done in a separate routine
        ! beforehand! 

        ! We then compute all intersections of this curve with all already
        ! available faces in the domain. Faces that are intersected in this way
        ! will be deleted and the new faces and intersections will be added.

        ! IMPORTANT
        !==========
        ! It should be noted that not all exceptional intersection cases are
        ! checked (see also the notes below for more information) so use this
        ! routine wisely. 

        ! Notes
        !======
        ! Note 1: if a curve part does not start or end in a new intersection or
        ! existing vertex, it will still be added as a face here. Clean-up should
        ! be done later on, as it is possible that after adding more curves,
        ! additional intersections are found with these segments, which can't be
        ! found if we remove them prematurely. 

        ! Note 2: we assume that if the start or end vertex of a curve/face is
        ! present, the first/last x, y coordinate of the curve is exactly equal to
        ! that vertex's coordinates. 

        ! Note 3: it is assumed that the curves form simple polygons

        ! Note 4: there may still be exceptional intersection cases which are not
        ! properly captured by this routine. To avoid these as much as possible, it
        ! is good practice to compute these special intersection points beforehand.
        ! Typically these are saddle points, extrema, tangency points, ... When
        ! starting the curve from this point, which is assumed here, only checks
        ! need to be done whether the starting or ending point coincides with an
        ! existing one. 

        ! Note 5: we now also return a logical vector that is true for
        ! faces that have been deleted in this routine from the original
        ! topomesh. Added faces are appended and can therefore also be 
        ! queried. 

        ! Algorithm
        !==========
        ! 
        ! 1) For each curve, we compute the intersections with all existing 
        !    topological faces. For each topological face, do:
        !       1.1) Compute intersections using standard polygon intersection routine
        !       1.2) All intersections that have been found are checked whether
        !       they are close (up to precision defined below) to an existing
        !       vertex. If they are not close, they are added as new vertices. If
        !       they are close, it is assumed that the intersection happens at one
        !       of the existing vertices and no new vertex is added. The
        !       intersection is then adjusted to be exactly this original vertex,
        !       and also the ID of the intersection is set to the ID of that
        !       vertex. 
        !       1.3) Now, we check the end points of the segments. If any end
        !       points coincide, we check the IDs of these points. If they are the
        !       same and non-zero, nothing must be done. If one is non-zero, then
        !       we update the zero ID to the non-zero one. If both are zero, we
        !       need to add a new vertex to the topology mesh. If they are
        !       non-zero, but not the same, we need to throw an error. 
        !       1.4) Add for each curve the start and end vertex also as
        !       intersection with updated ID. 
        ! 2) After computing all intersections and having partitioned the segments
        ! as stated above, we can remove all old faces and add the new ones. Note
        ! that, even if no intersections are found of an existing face, the
        ! algorithm above will add that face again as expected. This may result in
        ! some overhead, so it is best to call this function only once by
        ! precomputing all desired curves beforehand, if possible. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT), intent(in)      :: magneticField 
        type(ContourUDT), intent(inout)         :: contour
        integer(I8), intent(in)                 :: contourfsID, contourtype
        logical, allocatable, intent(out), optional     :: remfout(:)

        ! Auxiliary
        real(R8)                                :: dist, fsfval
        real(R8), allocatable, dimension(:)     :: xint, yint, s1r, s2r, &
            tfv, ts1r, tx, ty
        integer(I8)                             :: nint, fsID, nforig
        integer(I8), allocatable, dimension(:)  :: s1, s2, fID, vIDs, &
            vtypes, sortind, vf1, vf2, tvIDs, ts2, ts1
        logical                                 :: alreadyadded, &
            isinconsistent
        logical, allocatable, dimension(:)      :: iscoinciding, &
            delind, keepind
        type(PolygonUDT)                        :: cp 
        type(RealDynamicArrayUDT)               :: xda, yda, s1rda, &
            s2rda
        type(RealDynamicArrayUDT), allocatable  :: xfda(:), yfda(:)
        type(IntegerDynamicArrayUDT)            :: s1da, s2da, fda 

        ! Loop
        integer(I8)                             :: i, j, k, kold

        ! Initialize
        !===========
        ! Original number of faces
        nforig = topomesh%face%ntot

        ! Construct polygon from contour data
        call cp%Construct(contour%x, contour%y)

        ! Initialize dynamic arrays
        xda     = ConstructRealDynamicArray()
        yda     = ConstructRealDynamicArray()
        s1rda   = ConstructRealDynamicArray()
        s2rda   = ConstructRealDynamicArray()
        s1da    = ConstructIntegerDynamicArray()
        s2da    = ConstructIntegerDynamicArray()
        fda     = ConstructIntegerDynamicArray()

        ! Compute intersections 
        !======================
        ! Loop over all faces
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp shared(topomesh, fda, xda, yda, s1da, s2da, s1rda, s2rda, cp) & 
        !$omp private(i, xint, yint, s1, s2, s1r, s2r)
        do i  = 1, topomesh%face%ntot 
            ! Associate current face polygon 
            associate(&
                fp          => topomesh%face%pol(i),    &
                fpx         => topomesh%face%pol(i)%x,  &
                fpy         => topomesh%face%pol(i)%y)

            ! Compute intersections
            call PolygonIntersections(fp, cp, xint, yint, s1, s2, &
                s1r=s1r, s2r=s2r)

            ! Store
            !$omp critical
            call fda%Append(spread(i, 1, size(xint)))
            call xda%Append(xint)
            call yda%Append(yint)
            call s1da%Append(s1)
            call s2da%Append(s2)
            call s1rda%Append(s1r)
            call s2rda%Append(s2r)
            !$omp end critical

            ! Housekeeping
            end associate
        end do
        !$omp end parallel do

        ! Extract
        fID     = fda%Get()
        xint    = xda%Get()
        yint    = yda%Get()
        s1      = s1da%Get()
        s2      = s2da%Get()
        s1r     = s1rda%Get()
        s2r     = s2rda%Get()
        nint    = size(xint)

        ! Add vertices
        !=============
        ! Initialize the vertex IDs and vertex types
        allocate(vIDs(nint), vtypes(nint))
        vIDs = 0
        vtypes = 0

        ! Loop over all intersections
        do i = 1, nint 
            ! Unpack
            associate(&
                xinti       => xint(i),     &
                yinti       => yint(i),     &
                fIDi        => fID(i),      &
                s1ri        => s1r(i),      &
                s2ri        => s2r(i)       &
                )

            ! Initialize logicals
            alreadyadded    = .false. 
            isinconsistent  = .false. 

            ! Set vertex type (default: 0)
            if ((contourtype == TMfacebndID) .or. &
                (topomesh%face%type(fIDi) == TMfacebndID)) then 
                ! Boundary vertex but no tangency point
                vtypes(i) = TMvertexbndID
            end if 

            ! Check if this vertex is already added
            iscoinciding = (abs(xinti - topomesh%vert%x) < disttol) .and. &
                (abs(yinti - topomesh%vert%y) < disttol)
            if (any (iscoinciding) ) then 
                ! Sanity check
                if (count(iscoinciding) > 1) then 
                    ! This shouldn't happen if the vertices in topomesh
                    ! are unique. Throw error
                    call gdErrorHandler('InsertTopologicalMeshContour: ' // & 
                        'duplicate vertices seem to appear in topomesh, ' // &
                        'check input')
                end if 

                ! Vertex coincides up to disttol precision -  do 
                ! sanity checks and adjustments of vertex ID
                alreadyadded = .true. 

                ! Was the intersection in a start point, if yes, do 
                ! vertex IDs correspond? 
                if ((s1ri == 0_R8) .and. (s2ri == 0_R8)) then 

                    ! Intersection in start of face and start of contour
                    if (contour%startsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 1)) /= contour%startsaddle) .or. &
                            (topomesh%face%vert(fIDi, 1) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 1) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 1)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%startsaddle = topomesh%face%vert(fIDi, 1)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == 0_R8) .and. &
                    (s2ri == real(cp%ne, R8))) then 

                    ! Intersection in start of face and end of contour
                    if (contour%endsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 1)) /= contour%endsaddle) .or. &
                            (topomesh%face%vert(fIDi, 1) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 1) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 1)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%endsaddle = topomesh%face%vert(fIDi, 1)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%endsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == 0_R8)) then 

                    ! Intersection in end of face and start of contour    
                    if (contour%startsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 2)) /= contour%startsaddle) .or. &
                            (topomesh%face%vert(fIDi, 2) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 2) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 2)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%startsaddle = topomesh%face%vert(fIDi, 2)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%startsaddle

                elseif ((s1ri == real(topomesh%face%pol(fIDi)%ne, R8)) .and. &
                    (s2ri == real(cp%ne, R8))) then 

                    ! Intersection in end of face and end of contour
                    if (contour%endsaddle /=0) then 
                        if (((topomesh%face%vert(fIDi, 2)) /= contour%endsaddle) .or. &
                            (topomesh%face%vert(fIDi, 2) == 0)) then 
                            ! Inconsistent - throw error later
                            isinconsistent = .true.
                        end if
                    else 
                        if (topomesh%face%vert(fIDi, 2) /= 0) then 
                            ! Not necessarily an error, but do print a warning
                            ! as this isn't expected
                            print *, 'vertex ID: ', topomesh%face%vert(fIDi, 2)
                            print *, 'InsertTopologicalMeshContour: ' // & 
                                'contour saddle point corresponds to ' // & 
                                'topomesh vertex but was not identified as ' // &
                                'such. Adding ID and continuing...'
                            contour%endsaddle = topomesh%face%vert(fIDi, 2)
                        else
                            ! Inconsistent
                            isinconsistent = .true.
                        end if 
                    end if 

                    ! Set ID 
                    vIDs(i) = contour%endsaddle

                else 

                    ! Apparently we're lucky and we get an intersection
                    ! in a vertex that already exists but that is 
                    ! not part of the start or end of a face or 
                    ! contour. This is very unlikely so we throw a 
                    ! warning yet continue and add the vertex ID 
                    print *, 'InsertTopologicalMeshContour: ' // & 
                        'found intersection that exactly coincides with ' // & 
                        'existing mesh vertex, yet was not identified ' // & 
                        'previously (i.e. the intersection happens to be )' // & 
                        'haphazardly together with an existing mesh vertex). ' // & 
                        'proceeding to add vertex, but results may be suprising'

                    ! Set the vertex ID 
                    vIDs(i) = findloc(iscoinciding, .true., 1)

                end if

                ! Set type
                vtypes(i) = topomesh%vert%type(vIDs(i))

            end if
            
            ! Check if we encountered an inconsistency, if so - call 
            ! error and exit. 
            if (isinconsistent) then 
                call gdErrorHandler('InsertTopologicalMeshContour: ' // &
                    'inconsistency encountered when adding intersection ' // & 
                    'with face number: ', fIDi)
            end if 

            ! Add the vertex (if not already present)
            if (.not. alreadyadded) then 
                ! Check flux surface ID
                fsID = 0
                if (contourfsID /= 0 ) then 
                    ! Sanity check
                    if (topomesh%face%fsID(fID(i)) /= 0) then 
                        ! Normally contours shouldn't intersect... 
                        print *, 'InsertTopologicalMeshContour: contour with fsID ', &
                            contourfsID, ' intersects with face ', fID(i), &
                            ' with fsID ', topomesh%face%fsID(fID(i)), ' - unexpected.' // & 
                            ' Taking face ID'
                        
                        ! Take face ID for contour
                        fsID = topomesh%face%fsID(fID(i))
                    else 
                        ! Take contour ID
                        fsID = contourfsID
                    end if 
                else 
                    ! Take face ID
                    fsID = topomesh%face%fsID(fID(i))
                end if 

                allocate(tfv(1))
                call magneticField%interp%Evaluate([xinti], [yinti], 0, 0, tfv)
                call AddTopologicalMeshVertex(topomesh, xinti, yinti, &
                    tfv(1), vtypes(i), fsID)
                ! Add ID as well 
                vIDs(i) = topomesh%vert%ntot
                deallocate(tfv)
            end if 
                
            ! Housekeeping
            end associate
        end do

        ! Add contour faces
        !==================
        ! Sort intersections according to contour coordinate
        allocate(sortind(size(s1r)))
        call Sort(s2r, ind=sortind)
        s1r = s1r(sortind)
        fID = fID(sortind)
        s1 = s1(sortind)
        s2 = s2(sortind)
        xint = xint(sortind)
        yint = yint(sortind)
        vIDs = vIDs(sortind)
        deallocate(sortind)

        ! Hedge for duplicate intersections (possible with closed polygons 
        ! or if multiple faces intersect in the same point)
        allocate(keepind(size(s1r)))
        keepind = .true. 
        do j = 1, size(keepind)-1
            if ((s2r(j+1) - s2r(j)) == 0_R8) then 
                keepind(j+1) = .false. 
            end if 
        end do

        ! Add start and end points as intersections if they have an 
        ! ID (and if that ID is not already present as an intersection)
        tvIDs = pack(vIDs, keepind) 
        ts2 = pack(s2, keepind)
        deallocate(keepind)
        if ((contour%startsaddle /= 0)) then 
            if (size(vIDs) > 0) then ! need to hedge for zero intersections
                if (contour%startsaddle /= vIDs(1)) then 
                    tvIDs = [contour%startsaddle, tvIDs]
                    ts2 = [1, ts2]
                end if 
            else
                tvIDs = [contour%startsaddle, tvIDs]
                ts2 = [1, ts2]
            end if 
        end if 
        if ((contour%endsaddle /= 0)) then 
            if (size(vIDs) > 0) then ! need to hedge for zero intersections
                if (contour%endsaddle /= vIDs(size(vIDs))) then 
                    tvIDs = [tvIDs, contour%endsaddle]
                    ts2 = [ts2, cp%ne]
                end if 
            else
                tvIDs = [tvIDs, contour%endsaddle]
                ts2 = [ts2, cp%ne]
            end if 
        end if 

        ! Extract faces
        call ExtractTopologicalFacesFromPolygon(cp, tvIDs, ts2, topomesh%vert%x, &
            topomesh%vert%y, vf1, vf2, xfda, yfda)

        ! Hedge for too small faces
        allocate(keepind(size(xfda)))
        keepind = .true. 
        do i = 1, size(xfda)
            tx = xfda(i)%Get()
            ty = yfda(i)%Get()
            dist = sum(sqrt((tx(2:size(tx)) - tx(1:size(tx)-1))**2) + &
                (ty(2:size(ty)) - ty(1:size(ty)-1))**2)
            if (dist <= disttol) then 
                ! Remove
                print *, 'local face ID: ', i 
                print *, 'vertices: ', vf1(i), vf2(i)
                print *, 'InsertTopologicalMeshContour: not adding face ' // & 
                    'with vertex indices as mentioned above as it is ' // & 
                    'smaller than distance tolerance'
                keepind(i) = .false. 
            end if 
        end do 
        
        ! Add to faces
        do i = 1, size(xfda)
            if (keepind(i)) then 
                call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                    yfda(i), contourtype, contourfsID, contour%val)
            end if 
        end do 

        ! Housekeeping
        deallocate(keepind)
        
        ! Adjust existing faces
        !======================
        ! Sort intersections according to face index
        allocate(sortind(size(s1r)))
        call Sort(fID, ind=sortind)
        s2r = s2r(sortind)
        s1r = s1r(sortind)
        s1 = s1(sortind)
        s2 = s2(sortind)
        xint = xint(sortind)
        yint = yint(sortind)
        vIDs = vIDs(sortind)
        deallocate(sortind)

        ! Extract
        k = 0
        do while (k < size(fID)) 

            ! Update loop variables
            kold = k 
            k = findloc(fID, fID(kold+1), 1, back=.true.)

            ! Add start and end points as intersections if they have an 
            ! ID (and if that ID is not already present as an intersection)
            tvIDs = vIDs(kold+1:k) 
            ts1 = s1(kold+1:k)
            ts1r = s1r(kold+1:k)

            ! Sort along ts1r
            allocate(sortind(size(ts1r)))
            call Sort(ts1r, ind=sortind)
            tvIDs = tvIDs(sortind)
            ts1 = ts1(sortind)
            deallocate(sortind)

            ! Hedge for duplicate intersections
            allocate(keepind(size(tvIDs)))
            keepind = .true. 
            do j = 1, size(keepind)-1
                if ((ts1r(j+1) - ts1r(j)) == 0_R8)  then 
                    keepind(j+1) = .false. 
                end if 
                ! Also check vertex ID for non-closed polygons, since s1r may be off by machine precision 
                if ((tvIDs(j+1) == tvIDs(j)) .and. .not. topomesh%face%pol(fID(k))%isclosed) then 
                    keepind(j+1) = .false.
                end if 
            end do 
            tvIDs = pack(tvIDs, keepind) ! can simply reduce here, not used afterwards
            ts1 = pack(ts1, keepind)
            deallocate(keepind)

            if ((topomesh%face%vert(fID(k), 1) /= 0) .and. (topomesh%face%vert(fID(k), 1) /= tvIDs(1))) then 
                tvIDs = [topomesh%face%vert(fID(k), 1), tvIDs]
                ts1 = [1, ts1]
            end if 
            if ((topomesh%face%vert(fID(k), 2) /= 0) .and. (topomesh%face%vert(fID(k), 2) /= tvIDs(size(tvIDs)))) then 
                tvIDs = [tvIDs, topomesh%face%vert(fID(k), 2)]
                ts1 = [ts1, topomesh%face%pol(fID(k))%ne]
            end if 

            ! Extract faces
            call ExtractTopologicalFacesFromPolygon(&
                topomesh%face%pol(fID(k)), tvIDs, ts1, topomesh%vert%x, &
                topomesh%vert%y, vf1, vf2, xfda, yfda)

            ! Add to faces
            do i = 1, size(xfda)
                if (topomesh%face%fsID(fID(k)) /= 0) then 
                    fsfval = topomesh%fsfval%Get(topomesh%face%fsID(fID(k)))
                else
                    fsfval = 0.0_R8
                end if 
                call AddTopologicalMeshFace(topomesh, [vf1(i), vf2(i)], xfda(i), &
                    yfda(i), topomesh%face%type(fID(k)), topomesh%face%fsID(fID(k)), fsfval)

            end do 

        end do 

        ! Remove adjusted faces
        allocate(delind(topomesh%face%ntot))
        delind = .false. 
        delind(fID) = .true. 
        call RemoveTopologicalMeshFaceLogical(topomesh, delind)

        ! Return optional arguments
        if (present(remfout)) then 
            remfout = delind(1:nforig)
        end if 

    end subroutine 

    ! Topological mesh face splitting
    subroutine SplitTopologicalMeshFaces(topomesh, remfout)

        ! Description
        !============
        ! This routine checks the topological mesh for the following faces:
        !
        ! - faces with same vertex indices (assumed still different face)
        ! - faces with the same start and end vertex indices (closed faces)
        !
        ! These faces are split up in resp. two and three parts in order to arrive
        ! at a conventional mesh format, where each face has a unique set of vertex
        ! indices (regardless the order) and no faces that close upon themselves
        ! exist. This is a prerequisite when mesh cells are determined. 

        ! To split up the faces, we add vertex nodes with ID -1 at the splitting
        ! points. 

        ! Note: a logical is now returned indicating deleted faces in the
        ! original topomesh, if desired. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, allocatable, intent(out), optional     :: remfout(:)

        ! Auxiliary 
        integer(I8)                             :: nncf, nndf, np, &
            ind(1:2), nfinit, tf
        integer(I8), allocatable, dimension(:)  :: faceID
        real(R8)                                :: newvcfx(1:2), newvcfy(1:2), &
            newvcff(1:2), newvdfx(1), newvdfy(1), newvdff(1), fsfval
        
        logical, allocatable, dimension(:)      :: isclosedface, &
            isduplicateface, delind

        type(RealDynamicArrayUDT)               :: xrda, yrda 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Set initial amount of faces
        nfinit = topomesh%face%ntot 

        ! Unpack for ease (only for determination of which faces to delete)
        associate(face      => topomesh%face, &
            ntot            => topomesh%face%ntot)

        ! Check faces
        !============
        ! Closed faces
        isclosedface = face%vert(:, 1) == face%vert(:, 2)

        ! Open faces
        allocate(isduplicateface(ntot), faceID(ntot))
        isduplicateface = .false. 
        faceID = 0
        do i = 1, ntot-1
            do j = i+1, ntot
                ! Skip closed faces - should be dealt with separately, even if
                ! multiple are present
                if ((.not. isclosedface(i)) .and. (.not. isclosedface(j))) then 
                    ! Check vertices
                    if (any(face%vert(i, 1) == face%vert(j, :)) .and. &
                        any(face%vert(i, 2) == face%vert(j, :))) then 
                        isduplicateface(j) = .true.
                        faceID(j) = i
                    end if 
                end if
            end do 
        end do 

        ! Count
        nncf = count(isclosedface)
        nndf = count(isduplicateface)

        ! Issue messages
        if (any(isclosedface)) then 
            print *, 'SplitTopologicalMeshFaces: ', nncf, &
                ' closed faces detected, splitting up ...'
        end if 
        if (any(isduplicateface)) then 
            print *, 'SplitTopologicalMeshFaces: ', nndf, &
            ' faces with the same vertices detected, splitting up ...'
        end if 

        ! Determine new vertices
        !=======================
        ! Closed faces
        do i = 1, ntot
            if (isclosedface(i)) then 
                
                ! Get number of points of this face
                np = face%x(i)%Size()
                if (np < 4) then  ! end points should be the same and duplicate
                    ! Issue message: we cannot split up this boundary
                    call gdErrorHandler('SplitTopologicalMeshFaces: ' // & 
                        'closed face found with only two coordinates, ' // & 
                        'cannot split up')
                end if 
                
                ! Split up the face into parts with approx. equal number of
                ! vertices
                if (np == 4) then 
                    ind = [2, 3]
                else
                    ind = [floor(real(np, R8)/3_R8), ceiling(real(2*np, R8)/3_R8)]
                end if
                
                ! Get vertex coordinates
                newvcfx = face%x(i)%Get(ind)
                newvcfy = face%y(i)%Get(ind)

                ! Get flux surface value
                if (face%fsID(i) /= 0) then 
                    fsfval = topomesh%fsfval%Get(face%fsID(i))
                else
                    fsfval = 0.0_R8
                end if 
                newvcff = fsfval
                
                ! Insert new vertices
                do j = 1, 2
                    call AddTopologicalMeshVertex(topomesh, &
                        newvcfx(j), newvcfy(j), newvcff(j), &
                        TMvertexsplitID, face%fsID(i))
                end do       
                
                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot-1], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(1), ind(2))]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(1), ind(2))]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot-1, topomesh%vert%ntot], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

                ! Insert third face
                xrda = ConstructRealDynamicArray(face%x(i)%Get([(k, k = ind(2), face%x(i)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(i)%Get([(k, k = ind(2), face%y(i)%Size())]))
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(i), face%fsID(i), fsfval)

            elseif (isduplicateface(i)) then 

                ! Check which face to split - preference to aligned faces
                ! instead of boundary faces (these may be contour parts
                ! etc)

                if ((face%type(i) /= TMfacebndID) .and. (face%x(i)%Size() > 2)) then 
                    tf = i
                elseif ((face%type(faceID(i)) /= TMfacebndID) .and. (face%x(faceID(i))%Size() > 2)) then 
                    tf = faceID(i)
                elseif (face%x(i)%Size() > 2) then 
                    tf = i 
                else
                    tf = faceID(i)
                end if 

                ! Get number of points of this face
                np = face%x(tf)%Size()
                if (np < 3) then  ! end points should be the same and duplicate
                    ! Issue message: we cannot split up this boundary
                    print *, 'face vertices: ', face%vert(i, 1), face%vert(i, 2)
                    call WriteTopologicalMesh(topomesh,'topomesh_error')
                    call gdErrorHandler('SplitTopologicalMeshFaces: ' // & 
                        'face with same vertices found with only ' // & 
                        'two coordinates, cannot split up')
                end if 
                
                ! Split up the face into parts with approx. equal number of
                ! vertices. 
                ind(1) = np/2+1
                
                ! Get vertex coordinates
                newvdfx = face%x(tf)%Get(ind(1))
                newvdfy = face%y(tf)%Get(ind(1))

                ! Get flux surface value
                if (face%fsID(tf) /= 0) then 
                    fsfval = topomesh%fsfval%Get(face%fsID(tf))
                else
                    fsfval = 0.0_R8
                end if 
                newvcff = fsfval
                
                ! Insert new vertex
                call AddTopologicalMeshVertex(topomesh, newvdfx(1), &
                    newvdfy(1), newvdff(1), TMvertexsplitID, face%fsID(tf))

                ! Insert first face
                xrda = ConstructRealDynamicArray(face%x(tf)%Get([(k, k = 1, ind(1))]))
                yrda = ConstructRealDynamicArray(face%y(tf)%Get([(k, k = 1, ind(1))]))
                call AddTopologicalMeshFace(topomesh, &
                    [face%vert(i, 1), topomesh%vert%ntot], &
                    xrda, yrda, face%type(tf), face%fsID(tf), fsfval)

                ! Insert second face
                xrda = ConstructRealDynamicArray(face%x(tf)%Get([(k, k = ind(1), face%x(tf)%Size())]))
                yrda = ConstructRealDynamicArray(face%y(tf)%Get([(k, k = ind(1), face%y(tf)%Size())]))
                
                call AddTopologicalMeshFace(topomesh, &
                    [topomesh%vert%ntot, face%vert(i, 2)], &
                    xrda, yrda, face%type(tf), face%fsID(tf), fsfval)

            end if 

        end do

        ! Remove old faces
        !=================
        ! Set IDs correctly
        allocate(delind(face%ntot))
        delind = .false.
        delind(1:nfinit) = isclosedface .or. isduplicateface
        call RemoveTopologicalMeshFaceLogical(topomesh, delind)

        ! Set output
        if (present(remfout)) then 
            remfout = delind(1:nfinit)
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Topomesh trimming
    subroutine TrimTopologicalMesh(topomesh, magneticField, vessel, remfout)

        ! Description
        !============
        ! This function removes boundaries and vertices from a topological mesh if
        ! those vertices and boundaries lie outside of the domain. To determine the
        ! latter, the vessel levelset function is checked, which is
        ! positive when points are outside of the boundary. Any points that lie on
        ! the boundary should be identified as tangency points (type 4 or 5) or
        ! other type of boundary points (type 6).

        ! All other points are removed. We assume that all intersections of faces with
        ! the boundary have been computed and that all segments are added as
        ! separate faces.

        ! Faces that do not have a start and end vertex are removed per definition.
        ! If a start or end vertex has been deleted, also the faces that have this
        ! start or end vertex will be deleted. For other faces, it is checked if
        ! the majority of points lies inside the vessel, excluding start and end
        ! vertices (they may lie exactly on the vessel or should be deleted already
        ! if they lie outside of the vessel). If all vertices (except end points)
        ! lie outside the vessel, there is no ambiguity and the face is deleted.
        ! If the first node near the start/end points are still inside
        ! the vessel, but the rest outside, it is assumed that this is due to mesh
        ! accuracy (a message will be displayed). If more than one node on each
        ! side is detected, a warning is issued and the face is not deleted. This
        ! may be due to not having computed and inserted all intersections, or due
        ! to misuse of this routine. 

        ! Notes
        !======
        ! Note 1: actually, we could also check based on the start and end vertex
        ! type whether we should consider a certain boundary for deletion (actually
        ! only boundaries with 'regular' intersections or tangency points). Now, we
        ! simply check all boundaries. 

        ! Note 2: we rely on the face identifiers to not consider boundary faces.
        ! This is necessary, since boundary faces may lie just on or off the vessel
        ! contour, depending on how accurate it was traced. Therefore, this routine
        ! is likely to fail for these boundaries (as expected...). Boundary faces
        ! are expected to be of type 3. 

        ! Note 3: separatrix segments (type 4) that have both intersections
        ! in boundary faces are removed as well, since they shouldn't be
        ! critical for the topological mesh. 

        ! Note 4: we also return a logical index indicating which faces
        ! were removed in the original topological mesh if desired. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel 
        logical, allocatable, intent(out), optional     :: remfout(:)

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: Vv
        logical, allocatable, dimension(:)      :: outbnd, rmvert, &
            rmface

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            plf         => vessel%plfvessel,      &
            mfinterp    => magneticField%interp)

        ! Rebuild the vessel description to be sure
        !    allocate(bndpol(count((topomesh%face%type == TMfacebndID) .or. &
        !        (topomesh%face%type == TMfacealbndID))))
        !    bndpol = pack(topomesh%face%pol, (topomesh%face%type == TMfacebndID) .or. &
        !        (topomesh%face%type == TMfacealbndID))
        !    call bndps%Construct(bndpol)
        !    call ConstructVesselPolygonSet(vessel, bndps)
        !    allocate(PLF2DClosedExactOptionsUDT::bndplfoptions)
        !    call InitializePolygonLevelsetFunction2D(vessel%plfvessel, &
        !        vessel%polygonset, bndplfoptions)

        ! Vertices
        !=========
        ! Check if in boundary
        allocate(Vv(topomesh%vert%ntot))
        call plf%Evaluate(topomesh%vert%x, topomesh%vert%y, 0, 0, Vv)
        outbnd = Vv >= 0

        ! Check if we should remove it
        rmvert = outbnd .and. (topomesh%vert%type /= TMvertextp1ID) .and. &
            (topomesh%vert%type /= TMvertextp2ID) .and. (topomesh%vert%type /= TMvertexbndID)

        ! Remove these vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, rmvert)

        ! Faces
        !======
        ! Start by removing faces with zero start or end vertex
        rmface = (topomesh%face%vert(:, 1) == 0) .or. &
            (topomesh%face%vert(:, 2) == 0)
        
        do i = 1, topomesh%face%ntot
            ! Remove separatrix faces that have only intersections with
            ! boundary faces
            if (topomesh%face%type(i) == TMfacesepID) then 
                ! Set to true, will be set to false if saddle point present
                rmface(i) = .true. 
                if ((topomesh%face%vert(i, 1) /= 0) .and. (topomesh%face%vert(i, 2) /= 0)) then 
                    if ((topomesh%vert%type(topomesh%face%vert(i, 1)) /= TMvertexbndID) .or. &
                        (topomesh%vert%type(topomesh%face%vert(i, 2)) /= TMvertexbndID)) then 
                        rmface(i) = .false. 
                    end if 
                else
                    ! Make sure is removed because of zero vertex
                    rmface(i) = .true.
                end if 
            end if 

            ! Remove non-boundary faces that are outside of the domain
            if (all(topomesh%face%type(i) /= [TMfacebndID, TMfacealbndID])) then 
                ! Points in boundary (exclude end points)?
                deallocate(Vv)
                allocate(Vv(size(topomesh%face%x(i)%Get())))
                call plf%Evaluate(topomesh%face%x(i)%Get(), &
                    topomesh%face%y(i)%Get(), 0, 0, Vv)
                outbnd = Vv(2:size(Vv)-1) >= 0
                
                ! Check
                if (size(outbnd) == 0) then 
                    ! Face with only two vertices - only keep if both 
                    ! vertices are non-zero 
                    if (any(topomesh%face%vert(i, :) == 0)) then 
                        rmface(i) = .true. 
                    end if
                elseif (all(outbnd)) then 
                    ! Remove, no issue
                    rmface(i) = .true.
                elseif (all(outbnd(2:size(outbnd)-1)) .and. (size(outbnd) > 2)) then 
                    ! Remove, but display message
                    print *, 'face ID: ', i, 'face vertices: ', topomesh%face%vert(i, 1), topomesh%face%vert(i, 2)
                    print *, 'TrimTopologicalMesh: boundary removed ' // & 
                        'which still had first two points in domain'
                elseif ((.not. any(outbnd(2:size(outbnd)-1))) .and. (size(outbnd) > 2)) then 
                    ! Do not remove, nothing to do here
                else
                    ! Check if we simply had a very short boundary
                    if (size(outbnd) <= 2) then 
                        print *, 'TrimTopologicalMesh: face detected ' // & 
                            'with at most four points, may not be ' // & 
                            'properly diagnosed for removal. Try ' // & 
                            'increasing number of points by increasing ' //& 
                            'contour mesh resolution'
                        if (count(.not. outbnd) < count(outbnd)) then 
                            rmface(i) = .true.
                        end if
                    else
                        ! Points detected inside and outside - throw warning
                        print *, 'TrimTopologicalMesh: face detected ' // & 
                            'that has multiple points inside and outside ' // & 
                            'of domain. May result in unexpected outcome. Check input'

                        if (count(.not. outbnd) < count(outbnd)) then 
                            rmface(i) = .true.
                        end if 
                    end if 
                end if 
            end if 
        end do 

        ! Remove
        call RemoveTopologicalMeshFaceLogical(topomesh, rmface)

        ! Return optional output arguments
        if (present(remfout)) then 
            remfout = rmface
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Topomesh simplification
    subroutine SimplifyTopologicalMeshFaces(topomesh, remfout)

        ! Description
        !============
        ! This routine attempts to simplify the topological mesh by 
        ! merging faces (and hence deleting vertices) if they would
        ! form a single, unique face again. To achieve this, we loop 
        ! over all vertices and check the following conditions:
        ! - does the vertex only appear in two separate boundaries?
        ! - are those two boundaries of the same type?
        ! - is the vertex not a split vertex or extremum?
        ! If all these conditions are met, the neighbouring faces may
        ! be safely merged into a single face. Normally, this merging 
        ! shouldn't be necessary, unless e.g. separatrix parts are 
        ! removed during topological mesh trimming. This routine should
        ! therefore be called after adding all boundaries to the 
        ! topological mesh, but before adding cells and other data. 

        ! Note: we also allow deletion of tangency/boundary points IF
        ! they have only two neighbouring aligned vessel faces as boundaries. 
        ! If those boundaries would have different flux surface ID, 
        ! they are still merged. 

        ! Note: we now return a logical vector that is true for removed
        ! faces according to the old number of faces in the original 
        ! topomesh. New faces can then also be found, since they have
        ! the index count(remfout)+1:topomesh%face%ntot

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, allocatable, intent(out), optional     :: remfout(:)

        ! Auxiliary
        integer(I8)                             :: nforig, tf(1:2)
        integer(I8), allocatable, dimension(:)  :: rmvID, rmf1, rmf2, &
            fvert, oldfID
        logical, allocatable, dimension(:)      :: markv, markf, &
            appearstwice
        real(R8)                                :: fsfval
        real(R8), allocatable, dimension(:)     :: tempx, tempy 
        type(RealDynamicArrayUDT)               :: xda, yda 

        ! Loop 
        integer(I8)                             :: i, k 

        ! Initialize
        !===========
        ! Set the old face ID index
        oldfID = [(k, k = 1, topomesh%face%ntot)]

        ! Initialize the output if necessary
        if (present(remfout)) then 
            allocate(remfout(topomesh%face%ntot))
            remfout = .true. 
        end if 

        ! Keep looping until exit condition reached
        do while (.true.) 

            ! Store original sizes
            nforig = topomesh%face%ntot

            ! Mark vertices and faces for removal
            allocate(markv(topomesh%vert%ntot), markf(topomesh%face%ntot), &
                rmvID(0), rmf1(0), rmf2(0))
            markv = .false. 
            markf = .false.
            do i = 1, topomesh%vert%ntot
                ! Is it a tangency point, regular vertex, or 
                ! boundary vertex? 
                if (.not. any(topomesh%vert%type(i) == [TMvertexbndID, &
                    0, TMvertextp2ID, TMvertextp1ID])) then 
                    cycle
                end if 

                ! Does it only appear twice in face%vert?
                appearstwice = ((topomesh%face%vert(:, 1) == i) .or. (topomesh%face%vert(:, 2) == i))
                if (.not. (count(appearstwice) == 2)) then 
                    cycle 
                end if 
                tf = pack([(k, k = 1, topomesh%face%ntot)], appearstwice)

                ! Are both faces not yet marked for deletion?
                if (any(markf(tf))) then 
                    cycle ! do in a next iteration
                end if 

                ! Do both boundaries have the same type? 
                if (.not. (topomesh%face%type(tf(1)) == topomesh%face%type(tf(2)))) then 
                    cycle
                end if 

                ! Is the vertex a type 1 tangency point? Then we need to 
                ! check if both faces are aligned boundary faces (flux 
                ! surface ID doesn't matter)
                if (topomesh%vert%type(i) == TMvertextp1ID) then 
                    if (.not. ( (topomesh%face%type(tf(1)) == TMfacealbndID) &
                        .and. (topomesh%face%type(tf(2)) == TMfacealbndID))) then 
                        cycle 
                    end if 
                end if 

                ! Check if same flux surface ID - not necessary anymore? 
                !if (.not. (topomesh%face%fsID(tf(1)) == topomesh%face%fsID(tf(2)))) then 
                !    cycle
                !end if 

                ! If we got here, we passed all checks. Mark for merging
                ! and deletion
                markv(i) = .true.
                markf(tf) = .true. 
                rmvID = [rmvID, i]
                rmf1 = [rmf1, tf(1)]
                rmf2 = [rmf2, tf(2)]

            end do 

            ! Remove old vertex IDs
            oldfID = pack(oldfID, .not. markf(1:size(oldfID))) ! markf from size(oldfID)+1 are newly added faces of this routine

            ! Check exit condition
            if (count(markv) == 0) then 
                exit 
            end if 

            ! Merge faces
            !============
            do i = 1, size(rmf1)
                ! Check vertices
                if (topomesh%face%vert(rmf1(i), 2) == rmvID(i)) then 
                    ! First face is in good position 
                    xda = topomesh%face%x(rmf1(i))
                    yda = topomesh%face%y(rmf1(i))
                    if (topomesh%face%vert(rmf2(i), 1) == rmvID(i)) then 
                        ! Simply append 
                        fvert = [topomesh%face%vert(rmf1(i), 1), &
                            topomesh%face%vert(rmf2(i), 2)]
                        call xda%Append(topomesh%face%x(rmf2(i))%Get())
                        call yda%Append(topomesh%face%y(rmf2(i))%Get())

                    elseif (topomesh%face%vert(rmf2(i), 2) == rmvID(i)) then 
                        ! Need to flip second face
                        fvert = [topomesh%face%vert(rmf1(i), 1), &
                            topomesh%face%vert(rmf2(i), 1)]
                        tempx = topomesh%face%x(rmf2(i))%Get()
                        tempy = topomesh%face%y(rmf2(i))%Get()
                        call xda%Append(tempx(size(tempx):1:-1))
                        call yda%Append(tempy(size(tempy):1:-1))

                    else 
                        ! This is a bug
                        call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                            'bug detected')
                    end if 
                elseif (topomesh%face%vert(rmf1(i), 1) == rmvID(i)) then 
                    ! Need to flip first face
                    tempx = topomesh%face%x(rmf1(i))%Get()
                    tempy = topomesh%face%y(rmf1(i))%Get()
                    tempx = tempx(size(tempx):1:-1)
                    tempy = tempy(size(tempy):1:-1)
                    xda = ConstructRealDynamicArray(tempx)
                    yda = ConstructRealDynamicArray(tempy)
                    if (topomesh%face%vert(rmf2(i), 1) == rmvID(i)) then 
                        ! Simply append 
                        fvert = [topomesh%face%vert(rmf1(i), 2), &
                            topomesh%face%vert(rmf2(i), 2)]
                        call xda%Append(topomesh%face%x(rmf2(i))%Get())
                        call yda%Append(topomesh%face%y(rmf2(i))%Get())

                    elseif (topomesh%face%vert(rmf2(i), 2) == rmvID(i)) then 
                        ! Need to flip second face
                        fvert = [topomesh%face%vert(rmf1(i), 2), &
                            topomesh%face%vert(rmf2(i), 1)]
                        tempx = topomesh%face%x(rmf2(i))%Get()
                        tempy = topomesh%face%y(rmf2(i))%Get()
                        call xda%Append(tempx(size(tempx):1:-1))
                        call yda%Append(tempy(size(tempy):1:-1))
                        
                    else 
                        ! This is a bug
                        call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                            'bug detected')
                    end if 
                else 
                    ! This is a bug
                    call gdErrorHandler('SimplifyTopologicalMesh: ' // &
                        'bug detected')
                end if 

                ! Add the new face
                if (topomesh%face%fsID(rmf1(i)) /= 0) then 
                    fsfval = topomesh%fsfval%Get(topomesh%face%fsID(rmf1(i)))
                else
                    fsfval = 0.0_R8
                end if 
                call AddTopologicalMeshFace(topomesh, fvert, &
                    xda, yda, topomesh%face%type(rmf1(i)), &
                    topomesh%face%fsID(rmf1(i)), fsfval)

            end do 

            ! Extend the face deletion vertices
            markf = [markf, spread(.false., 1, topomesh%face%ntot - nforig)]

            ! Delete the faces
            call RemoveTopologicalMeshFaceLogical(topomesh, markf)

            ! Delete the vertices
            call RemoveTopologicalMeshVertexLogical(topomesh, markv)

            ! Recompute some required data
            call AddTopologicalMeshVertexFaces(topomesh)

            ! Data
            call AddTopologicalMeshData(topomesh)

            ! Housekeeping
            deallocate(markv, markf, rmvID, rmf1, rmf2)

        end do 

        ! Optional output arguments
        if (present(remfout)) then 
            remfout(oldfID) = .false. 
        end if 

    end subroutine 

    ! Garbage tangency point removal
    subroutine RemoveGarbageTangencyPoints(topomesh)

        ! Description
        !============
        ! This routine identifies 'garbage' tangency points and removes
        ! them by changing their type and applying the topomesh
        ! simplification algorithm (SimplifyTopologicalMeshFaces). 
        ! Garbage tangency points are defined as type 1 tangency points
        ! that are not limited on both sides by the same flux surface.
        ! This assumes that contours have been added to the topomesh!

        ! We now additionally remove boundary points that result in 
        ! garbage behavior related to tangency points. In particular, 
        ! we check type 2 tangency points and see if they have
        ! a face which is a boundary face that has a regular boundary
        ! point with the same flux surface ID. This should normally not
        ! happen and indicates that the tangency point contour originally
        ! intersected with the boundary, but that the intersection 
        ! removal was not successful. This happens very rarily though...

        ! Note: no additional interconnection data is updated

        ! Note: when aligned vessel parts are present, we can't check
        ! on flux surface ID alone but have to chain aligned faces from
        ! one side of the tangency point to (hopefully) the other. We do 
        ! this in a similar way as cells are formed (though of course we
        ! don't store all the data). If we end up from one radial face
        ! into the other, with only one aligned boundary in between, 
        ! then the tangency point is valid. Otherwise, if there are 
        ! multiple distinct aligned parts, the tangency point is removed. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh

        ! Auxiliary
        integer(I8)                             :: ntvfal, tv1, tv2, &
            nnonaligned, loc1, loc2
        integer(I8), allocatable, dimension(:)  :: tvf, &
            tvf1, tvf2, tv, tf
        logical, allocatable, dimension(:)      :: delf, delv

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Simplify to be sure
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Reconstruct vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Add data
        call AddTopologicalMeshData(topomesh)

        ! Recompute cells
        call AddTopologicalMeshCells(topomesh)

        ! Associate
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            cell        => topomesh%cell    &
            )

        ! Checks
        !=======
        ! Type 1 tangency point garbage
        !------------------------------
        ! Initialize
        allocate(delv(vert%ntot))
        delv = .false.

        ! Check which TPs should be removed
        do i = 1, cell%ntot
            ! Get all faces and vertices of this cell
            tv = cell%GetVert(i)
            tf = cell%GetFace(i)

            ! Is there a tangency point of type 1?
            if (all(vert%type(tv) /= TMvertextp1ID)) then 
                cycle
            end if 

            ! Are there more than two non-aligned faces?
            nnonaligned = 0
            do j = 1, size(TMfacenonalignedID)
                nnonaligned = nnonaligned + count(face%type(tf) == TMfacenonalignedID(j))
            end do 
            if (nnonaligned == 2) then 
                cycle
            end if 

            ! Are there less than two aligned faces?
            if (nnonaligned < 2) then 
                ! This is weird, print message but continue
                print *, 'RemoveGarbageTangencyPoints: cell detected with ' // &
                    'only one or no non-aligned faces, but with type 1 ' // & 
                    'tangency point. Unexpected, but continuing without ' // &
                    'removing tangency point'
                cycle 
            end if 

            ! At this point, there should be multiple aligned faces. If
            ! the tangency point is adjacent to two of them, and only 
            ! has those two faces as neighbours, remove it by retyping it
            ! as a regular boundary vertex
            do j = 1, size(tv)
                if (vert%type(tv(j)) == TMvertextp1ID) then 
                    ! Get faces
                    tvf = vert%GetFace(tv(j))

                    ! Check if it only has two faces
                    if (size(tvf) /= 2) then 
                        ! Skip -  we won't be able to remove this one
                        cycle 
                    end if 

                    ! Check if the two faces have the same type
                    if (face%type(tvf(1)) /= face%type(tvf(2))) then 
                        cycle 
                    end if 
                    
                    ! Check if the two faces are adjacent in the cell 
                    ! faces (the cell faces should be ordened)
                    loc1 = findloc(tf, tvf(1), 1, back=.false.)
                    loc2 = findloc(tf, tvf(2), 1, back=.false.)
                    if (loc1 == 0 .or. loc2 == 0) then 
                        ! Skip, something weird
                        cycle
                    elseif (loc1 - loc2 == 0) then 
                        ! Weird - should be same face then
                        print *, 'RemoveGarbageTangencyPoints: vertex ' // & 
                            'neighbouring faces appear to be the same. ' // & 
                            'Unexpected, but moving on...'
                        cycle
                    elseif ((abs(loc1 - loc2) == 1) .or. &
                        ((loc1 == 1) .and. (loc2 == size(tf))) .or. &
                        ((loc2 == 1) .and. (loc2 == size(tf)))) then 

                        ! Adjacent, remove
                        delv(tv(j)) = .true. 
                    else
                        cycle
                    end if 
                end if 
            end do 
        end do 

        ! 'Remove'
        where (delv) vert%type = TMvertexbndID

        ! Housekeeping
        deallocate(delv)

        ! Simplify 
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Reconstruct vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Type 2 tangency point garbage
        !------------------------------
        ! Initialize
        allocate(delv(vert%ntot), delf(face%ntot))
        delv = .false.
        delf = .false.

        ! Check which boundary tangency points should be removed
        do j = 1, vert%ntot
            if (vert%type(j) == TMvertextp2ID) then 
                ! Get the faces of this vertex
                tvf = vert%GetFace(j)

                ! Check if there are two faces with the same flux 
                ! surface ID. If not, we need to check boundary faces
                ntvfal = count(face%fsID(tvf) == vert%fsID(j))
                if (ntvfal < 2) then 
                    ! Get boundary faces, should be two
                    tvf = pack(tvf, face%type(tvf) == TMfacebndID)

                    ! Check
                    if (size(tvf) < 2) then 
                        ! This is possible and will normally not lead to
                        ! any bad behavior - skip
                        cycle
                    elseif (size(tvf) > 2) then 
                        ! This is not possible - throw error
                        call WriteTopologicalMesh(topomesh, 'topomesh_error')
                        print *, 'vertex: ', j 
                        call gdErrorHandler('RemoveGarbageTangencyPoints: ' // & 
                            'type 2 tangency point has more than two vessel ' // & 
                            'faces, unexpected')
                    end if 

                    ! Check the flux surface ID of non-tangency point vertices
                    ! of both faces
                    if (face%vert(tvf(1), 1) == j) then 
                        tv1 = face%vert(tvf(1), 2)
                    else
                        tv1 = face%vert(tvf(1), 1)
                    end if 
                    if (face%vert(tvf(2), 1) == j) then 
                        tv2 = face%vert(tvf(2), 2)
                    else
                        tv2 = face%vert(tvf(2), 1)
                    end if 

                    ! Check first vertex
                    if (vert%fsID(tv1) == vert%fsID(j)) then 
                        ! Mark for deletion
                        delv(tv1) = .true.
                        delf(tvf(1)) = .true.

                        ! Print
                        print *, 'RemoveGarbageTangencyPoints: removing ', &
                            'type 2 tangency point vertex ', tv1

                        ! Get other faces of this vertex
                        tvf1 = vert%GetFace(tv1)
                        tvf1 = pack(tvf1, tvf1 /= tvf(1))

                        ! Adjust the vertex of these faces
                        do k = 1, size(tvf1)
                            ! Adjust the end point
                            if (face%vert(tvf1(k), 1) == tv1) then 
                                face%vert(tvf1(k), 1) = j 
                                call face%x(tvf1(k))%Set(1, vert%x(j))
                                call face%y(tvf1(k))%Set(1, vert%y(j))
                            else
                                face%vert(tvf1(k), 2) = j 
                                call face%x(tvf1(k))%Set(face%x(tvf1(k))%Size(), vert%x(j))
                                call face%y(tvf1(k))%Set(face%y(tvf1(k))%Size(), vert%y(j))
                            end if

                            ! Reconstruct the polygon
                            call face%pol(tvf1(k))%Construct(face%x(tvf1(k))%Get(), &
                                    face%y(tvf1(k))%Get())
                        end do 
                    end if

                    ! Check second vertex
                    if (vert%fsID(tv2) == vert%fsID(j)) then 
                        ! Mark for deletion
                        delv(tv2) = .true.
                        delf(tvf(2)) = .true.

                        ! Print
                        print *, 'RemoveGarbageTangencyPoints: removing ', &
                            'vertex ', tv2

                        ! Get other faces of this vertex
                        tvf2 = vert%GetFace(tv2)
                        tvf2 = pack(tvf2, tvf2 /= tvf(2))

                        ! Adjust the vertex of these faces
                        do k = 1, size(tvf2)
                            ! Adjust the end point
                            if (face%vert(tvf2(k), 1) == tv2) then 
                                face%vert(tvf2(k), 1) = j 
                                call face%x(tvf2(k))%Set(1, vert%x(j))
                                call face%y(tvf2(k))%Set(1, vert%y(j))
                            else
                                face%vert(tvf2(k), 2) = j 
                                call face%x(tvf2(k))%Set(face%x(tvf2(k))%Size(), vert%x(j))
                                call face%y(tvf2(k))%Set(face%y(tvf2(k))%Size(), vert%y(j))
                            end if

                            ! Reconstruct the polygon
                            call face%pol(tvf2(k))%Construct(face%x(tvf2(k))%Get(), &
                                    face%y(tvf2(k))%Get())
                        end do 
                    end if
                end if 
            end if 
        end do

        ! Remove
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)


        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Boundary split vertex contours
    subroutine AddBoundarySplitVertexContours(topomesh, magneticField, vessel, &
        fieldtracer)

        ! Description
        !============
        ! This routine adds contours that start from a split vertex 
        ! located on a boundary. These are necessary to deal with 
        ! boundaries that have been split and would result in 
        ! non-conforming flux tubes (perhaps expanding the flux tube
        ! concept will be necessary at some point...)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(MagneticFieldUDT),intent(in)       :: magneticField 
        type(VesselUDT), intent(inout)          :: vessel
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 

        ! Auxiliary
        real(R8)                                :: dist
        integer(I8), allocatable, dimension(:)  :: allcurvetypes, allfsIDs
        logical, allocatable, dimension(:)      :: keepind

        type(ContourUDT), allocatable           :: tc(:), allc(:), &
            allbsvc(:)
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        type(IntegerDynamicArrayUDT)            :: curvetypes, fsIDs


        ! Loop 
        integer(I8)                             :: i, k
        
        ! Initialize
        !===========
        ! Associate
        associate(nfs       => topomesh%nfs)

        ! Initialize the contour structure & dynamic arrays
        allocate(allc(0))
        curvetypes = ConstructIntegerDynamicArray()
        fsIDs = ConstructIntegerDynamicArray()

        ! First, trace all contours
        allocate(allbsvc(0))
        !$omp parallel do default(none) schedule(dynamic) &
        !$omp private(i, k, tc, keepind, dist) & 
        !$omp shared(topomesh, fieldtracer, allbsvc, curvetypes, &
        !$omp fsIDs) 
        do i = 1, topomesh%vert%ntot
            if ((topomesh%vert%type(i) == TMvertexsplitID) .and. &
                topomesh%vert%fsID(i) == 0) then 
                ! Trace contour
                tc = fieldtracer%TraceContours([topomesh%vert%x(i)], [topomesh%vert%y(i)])

                ! Process
                call CleanContours(tc)

                ! Checks
                allocate(keepind(size(tc)))
                keepind = .true.
                do k = 1, size(tc)
                    ! Is the starting point the actual given start point? (should
                    ! be exactly the same since added as starting point in the
                    ! contouring algorithm)
                    dist = sqrt((tc(k)%x(1) - topomesh%vert%x(i))**2 + &
                        (tc(k)%y(1) - topomesh%vert%y(i))**2)
                    if (dist > 0.0_R8) then 
                        print *, 'AddTopologicalMEshContours: split vertex ' // & 
                            'contour segment found that does not start ' // & 
                            'in given tangency point. Removing...'
                        keepind(k) = .false.
                    end if
                    if (tc(k)%startsaddle /= i) then
                        ! Normally this should come from the contour 
                        ! tracer, but we can add it afterwards as well 
                        print *, 'AddTopologicalMeshContours: split vertex ' // &
                            'contour segment found that starts in given ' // & 
                            'point, but that does not have the ' // &
                            'starting saddle point ID as tangency point. ' // & 
                            'adjusting starting ID...'
                        tc(k)%startsaddle = i
                    end if 
                    if (tc(k)%isclosed .and. (tc(k)%endsaddle /= i)) then 
                        ! Ensure start and end saddle point are the same
                        tc(k)%endsaddle = i
                    end if 
                end do

                ! Remove
                tc = pack(tc, keepind)
                deallocate(keepind)

                ! Add
                !$omp critical
                allbsvc = [allbsvc, tc]
                call curvetypes%Append(spread(TMfacepolID, 1, size(tc)))
                
                ! Add flux surface ID
                topomesh%nfs = topomesh%nfs + 1
                topomesh%vert%fsID(i) = topomesh%nfs
                topomesh%vert%type(i) = TMvertexbndID
                topomesh%vert%fval(i:i) = fieldtracer%Evaluate([topomesh%vert%x(i)], [topomesh%vert%y(i)])
                call fsIDs%Append(spread(topomesh%nfs, 1, size(tc)))
                !$omp end critical
            end if 
        end do 
        !$omp end parallel do 

        ! Add contours (IDs etc should be done already...)
        allc = [allc, allbsvc]

        ! Add to the topology mesh
        !=========================
        ! Clean the contours (just to be sure)
        if (size(allc) > 0) then 
            call CleanContours(allc)

            allocate(allpc(size(allc)))
            do i = 1, size(allc)
                call allpc(i)%Construct(allc(i)%x, allc(i)%y)
            end do 
            call tempps%Construct(allpc)
            call tempps%WriteData('topocontours_all_beforeinsertion2')

            ! Add the contours
            allcurvetypes = curvetypes%Get()
            allfsIDs = fsIDs%Get()
            do i = 1, size(allc)
                call InsertTopologicalMeshContour(topomesh, magneticField, &
                    allc(i), allcurvetypes(i), allfsIDs(i))
            end do 

            ! Trim the topological mesh
            call TrimTopologicalMesh(topomesh, magneticField, vessel)

            ! Split boundaries
            call SplitTopologicalMeshFaces(topomesh)   
        end if
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Topomesh cleaner
    subroutine CleanTopologicalMesh(topomesh)

        ! Description
        !============
        ! This routine cleans up the topological mesh without further 
        ! modifying it (this routine should only be called at the 
        ! end of topomesh construction). At the moment, this only 
        ! includes removal of flux surface IDs that are not present 
        ! anymore in any of the faces or vertices. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)           :: topomesh

        ! Auxiliary
        logical, allocatable, dimension(:)          :: keepind
        integer(I8)                                 :: nfsold
        integer(I8), allocatable, dimension(:)      :: fsIDold, fsIDmap
        real(R8), allocatable, dimension(:)         :: fsfvalold

        ! Loop
        integer(I8)                                 :: i, k 

        ! Initialize
        !===========
        ! Unpack
        nfsold      = topomesh%nFs
        fsIDold     = topomesh%fsID%Get()
        fsfvalold   = topomesh%fsfval%Get()

        ! Initialize
        allocate(fsIDmap(nfsold), keepind(nfsold))
        keepind = .false. ! will be set to true if found
        fsIDmap = 0_I8

        ! Determine used flux surface IDs
        !================================
        ! Used by vertices
        do i = 1, topomesh%vert%ntot
            if (topomesh%vert%fsID(i) /= 0) then 
                keepind(topomesh%vert%fsID(i)) = .true.
            end if 
        end do 

        ! Used by faces
        do i = 1, topomesh%face%ntot
            if (topomesh%face%fsID(i) /= 0) then 
                keepind(topomesh%face%fsID(i)) = .true.
            end if 
        end do 

        ! Construct mapping
        !==================
        k = 0
        do i = 1, nfsold
            if (keepind(i)) then 
                k = k + 1
                fsIDmap(i) = k
            else
                fsIDmap(i) = 0
            end if
        end do 

        ! Remap
        !======
        ! Vertices
        do i = 1, topomesh%vert%ntot
            if (topomesh%vert%fsID(i) /= 0) then 
                topomesh%vert%fsID(i) = fsIDmap(topomesh%vert%fsID(i))
            end if 
        end do 

        ! Faces
        do i = 1, topomesh%face%ntot
            if (topomesh%face%fsID(i) /= 0) then 
                topomesh%face%fsID(i) = fsIDmap(topomesh%face%fsID(i))
            end if 
        end do 

        ! Delete
        topomesh%nFs = maxval(fsIDmap)
        topomesh%fsID = ConstructIntegerDynamicArray([(k, k = 1, topomesh%nFs)])
        topomesh%fsfval = ConstructRealDynamicArray(pack(topomesh%fsfval%Get(), keepind))

    end subroutine

    !------------------------------------------------------------------!
    !                  BASIC TOPOLOGICAL MESH ADAPTATIONS              !
    !------------------------------------------------------------------!
    ! Core boundary contours
    subroutine AddTopologicalMeshCoreBoundaries(topomesh, magneticField, &
        vessel, fieldtracer, options)

        ! Description
        !============
        ! This routine adds additional poloidal contours (type 2) for core parts.
        ! True limiter configurations are not yet supported (a warning is thrown). 
        ! The field value of this contour is determined as frac*(psiO - psiX) +
        ! psiX, where psiX is the value of the field at the X-point (or other 
        ! point), and psiO the value at the extremum (i.e. frac = 0 -> x-point,
        ! frac = 1 -> extremum). Note that values of frac close to 0 or 1 may lead
        ! to problems further downstream, as these contours may be very coarse
        ! depending on the tracing grid size. 

        ! It is assumed that all other necessary contours have been added (i.e.
        ! this routine should be invoked after 'AddTopologicalMeshContours).
        ! Otherwise, the introduced core boundary may not be properly determined.
        ! The remainder of this routine is much alike AddTopologicalMeshContours.

        ! Note: the type of the flux surface is currently a regular
        ! poloidal face type (i.e. TMfacepolID) instead of a core type.
        ! This is because core boundaries are only actual outer 
        ! boundaries until the core regions have been removed. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(inout)          :: vessel 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, s1, s2, &
            contourtypes, fsIDs
        real(R8)                                :: thisfval, traceval 
        real(R8), allocatable, dimension(:)     :: xint, yint 
        logical, allocatable, dimension(:)      :: keepind 
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            vert            => topomesh%vert    &
        )

        ! Allocate
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! Trace contours
        !===============
        ! Loop
        do i = 1, vert%ntot 
            ! Check if we should trace for this extremum
            if (.not. ((vert%type(i) == TMvertexminID) .or. (vert%type(i) == TMvertexmaxID))) then 
                cycle 
            end if 

            ! Get the faces that have this extremum
            tf = pack([(k, k = 1, topomesh%face%ntot)], any(topomesh%face%vert == i, dim=2))
            !tf = findloc(any(topomesh%face%vert == i, dim=2), .true.)

            ! Sanity check
            if (size(tf) == 0) then 
                ! No faces found - limiter-like configuration, not yet 
                ! supported here
                print *, 'AddTopologicalMeshCoreBoundaries: limiter-like case ' // & 
                    'detected which is not yet supported. Not adding core ' // & 
                    'boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
                cycle 
            end if 

            ! Choose face with closest field value
            if (allocated(tfv)) then 
                deallocate(tfv)
            end if 
            allocate(tfv(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then 
                    tfv(j) = topomesh%face%vert(tf(j), 2)
                else
                    tfv(j) = topomesh%face%vert(tf(j), 1) 
                end if 
            end do 
            thisf = tf(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1))
            thisfval = vert%fval(tfv(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1)))

            ! Compute the field value to trace
            traceval = options%coreboundariesfrac*(vert%fval(i) - thisfval) + thisfval

            ! Trace contours 
            tc = fieldtracer%TraceContours([traceval])
            call CleanContours(tc)

            ! Check which contour is closed and intersects with the 
            ! current face
            allocate(keepind(size(tc)))
            keepind = .true. 
            do j = 1, size(tc)
                ! Convert to polygon
                call tcp%Construct(tc(j)%x, tc(j)%y)

                ! Compute intersections
                call PolygonIntersections(topomesh%face%pol(thisf), tcp, &
                    xint, yint, s1, s2)

                ! Check
                if (size(xint) == 0) then 
                    keepind(j) = .false. 
                end if 
            end do

            ! Remove contours without intersections
            tc = pack(tc, keepind)
            deallocate(keepind)

            ! Check
            if (size(tc) == 0) then 
                ! No boundaries found
                print *, 'AddTopologicalMeshCoreBoundaries: no ' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'Not adding any core boundaries for extremum with ' // & 
                    'coordinates: (', vert%x(i), ', ', vert%y(i), ')'
                cycle 
            elseif (size(tc) > 1) then 
                ! Multiple boundaries found - also not expected, but
                ! may not be a problem 
                print *, 'AddTopologicalMeshCoreBoundaries: multiple' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'May lead to too many core boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
            end if 

            ! Concatenate
            allc = [allc, tc]
            contourtypes = [contourtypes, spread(TMfacepolID, 1, size(tc))]

            ! Add flux surface ID
            nfs = nfs + 1
            fsIDs = [fsIDs, spread(nfs, 1, size(tc))]

        end do 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine 

    ! PF boundary contours
    subroutine AddTopologicalMeshPFBoundaries(topomesh, magneticField, &
        vessel, fieldtracer, options)

        ! Description
        !============
        ! This routine adds 'private flux' boundaries where applicable
        ! (here, the PF is defined as a region that has as a tangency 
        ! point with contour lying outside of the vessel, i.e. type 1 
        ! tangency points, and of which at least one of the two 
        ! neighbouring faces has a vertex that belongs to a separatrix
        ! segment). 
         
        ! The field value of this contour is determined as frac*(psiTP - psiX) +
        ! psiX, where psiX is the value of the field at the X-point (or other 
        ! point), and psiO the value at the tangency point (i.e. frac = 0 -> x-point,
        ! frac = 1 -> tp). Note that values of frac close to 0 or 1 may lead
        ! to problems further downstream, as these contours may be very coarse
        ! depending on the tracing grid size. 

        ! It is assumed that all other necessary contours have been added (i.e.
        ! this routine should be invoked after 'AddTopologicalMeshContours).
        ! Otherwise, the introduced boundary may not be properly determined.
        ! The remainder of this routine is much alike AddTopologicalMeshContours.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(inout)          :: vessel 
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, s1, s2, &
            contourtypes, fsIDs, sepfsIDs
        real(R8)                                :: thisfval, traceval
        real(R8), allocatable, dimension(:)     :: xint, yint
        logical                                 :: skipvertex
        logical, allocatable, dimension(:)      :: keepind, hasfv
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 
        ! real(R8), allocatable, dimension(:)     :: 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            vert            => topomesh%vert    &
        )

        ! Allocate
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! Get all separatrix flux surface IDs
        sepfsIDs = topomesh%GetSeparatrixFluxSurfaceIDs()

        ! Trace contours
        !===============
        ! Loop
        do i = 1, vert%ntot 
            ! Check if we should trace for this point
            if (.not. ((vert%type(i) == TMvertextp1ID))) then 
                cycle 
            end if 

            ! Get the faces that have this extremum
            hasfv = any(topomesh%face%vert == i, dim=2)
            if (allocated(tf)) then 
                deallocate(tf)
            end if 
            allocate(tf(count(hasfv)))
            tf = pack([(k, k = 1, topomesh%face%ntot)], hasfv)

            ! Sanity check
            if (size(tf) /= 2) then 
                ! Weird - normally, each tangency point of this type
                ! should only have two (radial) faces - ignore
                print *, 'AddTopologicalMeshPFBoundaries: tangeny point ' // & 
                    'vertex number: ', i, 'does not have exactly two ' // & 
                    'faces - unexpected. Skipping this tangency point...'
                cycle 
            end if 

            ! Check if both faces are target faces
            if (allocated(tfv)) then 
                deallocate(tfv)
            end if 
            allocate(tfv(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then 
                    tfv(j) = topomesh%face%vert(tf(j), 2)
                else
                    tfv(j) = topomesh%face%vert(tf(j), 1) 
                end if 
            end do 

            skipvertex = .false.
            do j = 1, size(tfv)
                if (.not. any(topomesh%vert%fsID(tfv(j)) == sepfsIDs)) then 
                    skipvertex = .true. 
                    exit
                end if 
            end do 
            if (skipvertex) then 
                cycle
            end if 

            ! If we got here, check which value comes closest to tp value
            thisf = tf(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1))
            thisfval = vert%fval(tfv(minloc(abs(vert%fval(i) - vert%fval(tfv)), dim=1)))

            ! Compute the field value to trace
            traceval = options%PFboundariesfrac*(vert%fval(i) - thisfval) + thisfval

            ! Trace contours 
            tc = fieldtracer%TraceContours([traceval])
            call CleanContours(tc)

            ! Check which contour intersects with both boundaries
            allocate(keepind(size(tc)))
            keepind = .true. 
            do j = 1, size(tc)
                ! Convert to polygon
                call tcp%Construct(tc(j)%x, tc(j)%y)

                ! Compute intersections with faces
                do k = 1, size(tf)
                    call PolygonIntersections(topomesh%face%pol(tf(1)), tcp, &
                        xint, yint, s1, s2)

                    ! Check
                    if (size(xint) == 0) then 
                        keepind(j) = .false. 
                    end if 
                end do 
            end do

            ! Remove contours without intersections
            tc = pack(tc, keepind)
            deallocate(keepind)

            ! Check
            if (size(tc) == 0) then 
                ! No boundaries found
                print *, 'AddTopologicalMeshPFBoundaries: no ' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'Not adding any core boundaries for extremum with ' // & 
                    'coordinates: (', vert%x(i), ', ', vert%y(i), ')'
                cycle 
            elseif (size(tc) > 1) then 
                ! Multiple boundaries found - also not expected, but
                ! may not be a problem 
                print *, 'AddTopologicalMeshPFBoundaries: multiple' // & 
                    'contours found that intersect with the extremum face. ' // & 
                    'May lead to too many core boundaries for extremum with coordinates (', & 
                    vert%x(i), ' , ', vert%y(i), ')'
            end if 

            ! Concatenate
            allc = [allc, tc]
            contourtypes = [contourtypes, spread(TMfacepolID, 1, size(tc))]

            ! Add flux surface ID
            nfs = nfs + 1
            fsIDs = [fsIDs, spread(nfs, 1, size(tc))]

            ! Housekeeping
            deallocate(tf, tfv)

        end do 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine 

    ! Limiter region removal
    subroutine RemoveTopologicalMeshLimiterRegions(topomesh, &
        magneticField, fieldtracer, streamlinetracer, options)

        ! Description
        !============
        ! This routine removes limiter-like regions (i.e. regions that 
        ! form disc cells, but that do not have any radial face, so 
        ! basically a minimum/maximum without radial line) by 
        ! introducing a connection between an extremum and a limiting
        ! tangency point. This connection is a streamline if it is
        ! found, otherwise it is a straight line between extremum and
        ! tangency point (this may perform poorly though, a warning will
        ! be issued). It is assumed that the necessary contours have 
        ! already been added, but that no core or other unnecessary 
        ! contours are present, nor that any regions have been removed. 
        ! We do assume that all remaining points are within the vessel
        ! boundary. 

        ! The algorithm is as follows: 
        ! 1) check which vertices are extrema without any faces 
        ! -> these are extrema in limiter-like configurations. 
        ! 2) retrace from these extrema a df/dx line (or df/dy line) 
        ! 3) check with which tangency point contour this line crosses 
        ! first -> that is the limiting tangency point
        ! 4) trace from the tangency point a streamline to the extremum
        ! 5) add that line as a face to the topological mesh

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        type(MagneticFieldUDT), intent(in)      :: magneticField
        class(ContourTracerUDT), intent(inout)  :: fieldtracer
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer
        class(TopomeshOptionsUDT), intent(in)   :: options 

        ! Auxiliary
        logical                                 :: hasfaces
        logical, allocatable, dimension(:)      :: islimiterextremum
        integer(I8)                             :: ngp, nle, faceID, &
            thistp
        integer(I8), allocatable, dimension(:)  :: sface, sfc, tfID, &
            tpfsID
        real(R8)                                :: xb(1:2), yb(1:2)
        real(R8), allocatable, dimension(:)     :: fx, fy, xg, yg, &
            xint, yint, srface, srfc, tsrfc
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer
        type(ContourUDT), allocatable           :: fxc(:), fyc(:), fc(:)
        type(StreamlineUDT), allocatable        :: streamlines(:)
        type(IntegerDynamicArrayUDT)            :: intfID
        type(RealDynamicArrayUDT)               :: intsrfc
        class(RealDynamicArrayUDT), allocatable :: tx, ty
        class(StreamlineTracerUDT), allocatable :: mystreamlinetracer

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            vert        => topomesh%vert, &
            face        => topomesh%face)        

        ! Initialize
        allocate(islimiterextremum(vert%ntot))
        islimiterextremum = .false. 

        ! Tracer bounds
        xb(1:2) = [minval(magneticField%R), maxval(magneticField%R)]
        yb(1:2) = [minval(magneticField%Z), maxval(magneticField%Z)]

        ! Copy and adjust tracer
        mystreamlinetracer = streamlinetracer 

        ! Set step size
        mystreamlinetracer%step = 0.1
        mystreamlinetracer%nsteps = 2000

        ! Determine limiter exterma
        !==========================
        ! Loop over all vertices
        do i = 1, vert%ntot 
            ! Check if extremum
            if (.not. any(topomesh%vert%type(i) == [TMvertexmaxID, TMvertexminID])) then 
                ! Skip
                cycle 
            end if 

            ! Check if this extremum has faces
            j = 1
            hasfaces = .false. 
            do while (j <= face%ntot)
                if (any(face%vert(j, :) == i)) then 
                    hasfaces = .true.
                    exit
                else
                    j = j + 1
                end if 
            end do 
            if (hasfaces) then 
                ! Skip
                cycle 
            end if 

            ! Extremum without faces found, set to true
            islimiterextremum(i) = .true.
        end do 

        ! Check if we should continue
        nle = count(islimiterextremum)
        if (nle == 0) then 
            return 
        end if 

        ! Initialize tracers
        !===================
        ! If we got here, limiter extrema were found. Initializing tracers
        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Evaluate
        allocate(fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers - set only limiter extrema as saddle points
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)

        ! Construct faces
        !================
        ! Get all flux surface IDs that belong to type 1 tangency points
        allocate(tpfsID(count(vert%type == TMvertextp2ID)))
        tpfsID = pack(vert%fsID, vert%type == TMvertextp2ID)

        ! Loop over all vertices
        do i = 1, vert%ntot
            ! Check if this is a limiter extremum
            if (.not. islimiterextremum(i)) then 
                ! Skip
                cycle
            end if 

            ! Trace the contour lines from this point
            fxc = fxtracer%TraceContours([vert%x(i)], [vert%y(i)]) ! not exactly df/dx = 0 but ok
            fyc = fytracer%TraceContours([vert%x(i)], [vert%y(i)])

            ! Hedge for possible non-allocation 
            if (.not. allocated(fxc)) then 
                allocate(fxc(0))
            end if 
            if (.not. allocated(fyc)) then 
                allocate(fyc(0))
            end if

            ! Keep only lines that start/end in this 

            ! Check
            if (size(fxc) == 0 .and. size(fyc) == 0) then 
                ! Strange, no lines found
                call gdErrorHandler('RemoveTopologicalMeshLimiterRegions: ' // &
                    'could not find any df/dx=0 or df/dy=0 lines in extremum, ' // &
                    'unexpected. Could be a bug')
            end if 
            fc = [fxc, fyc]

            ! For the first of these lines, find intersections with any 
            ! tangency point contour lines
            intfID = ConstructIntegerDynamicArray()
            intsrfc = ConstructRealDynamicArray()
            do j = 1, face%ntot
                if (any(face%fsID(j) == tpfsID)) then 
                    ! Look for intersections
                    call SimplePolygonIntersections(face%x(j)%Get(), &
                        face%y(j)%Get(), fc(1)%x, fc(1)%y, xint, yint, &
                        sface, sfc, srface, srfc)

                    ! Check if any found
                    if (size(xint) > 0) then 
                        ! Append 
                        call intfID%Append(j)
                        call intsrfc%Append(minval(srfc)) ! just take minimal value
                    end if 
                end if 
            end do 

            ! If no intersections found -> throw error
            if (intfID%Size() == 0) then 
                call gdErrorHandler('RemoveTopologicalMeshLimiterRegions: ' // & 
                    'no intersections with any tangency point contour found, this ' // & 
                    'may be a bug')
            end if 

            ! Find intersection closest to extremum (i.e. the one with 
            ! smallest segment index)
            tsrfc = intsrfc%Get()
            tfID = intfID%Get()
            faceID = tfID(minloc(tsrfc, 1))

            ! Find the corresponding tangency point
            thistp = findloc(vert%type == TMvertextp2ID .and. &
                vert%fsID == face%fsID(faceID), .true., 1)

            ! Trace a contour line from this point to the extremum
            if (vert%type(i) == TMvertexminID) then 
                ! Need to go downhill, so against the gradient
                streamlines = mystreamlinetracer%TraceStreamlines([vert%x(thistp)], &
                    [vert%y(thistp)], xb, yb, [-1_I8])
            else
                ! Need to go uphill, with the gradient
                streamlines = mystreamlinetracer%TraceStreamlines([vert%x(thistp)], &
                    [vert%y(thistp)], xb, yb, [1_I8])
            end if 

            ! Checks
            if (size(streamlines) /= 1) then 
                ! Print warning and simply connect tangency point and extremum
                print *, 'extremum vertex: ', i, 'tangeny vertex: ', thistp
                print *, 'RemoveTopologicalMeshLimiterRegions: ' // & 
                    'could not trace streamline from tangency point to ' // &
                    'extremum, adding straight line between these points'

                ! Add straight line
                tx = ConstructRealDynamicArray([vert%x(thistp), vert%x(i)])
                ty = ConstructRealDynamicArray([vert%y(thistp), vert%y(i)])
                call AddTopologicalMeshFace(topomesh, [thistp, i], &
                    tx, ty, TMfaceradID, 0_I8, 0.0_R8)
            else
                ! Append the extremum to the streamline
                tx = ConstructRealDynamicArray([streamlines(1)%x, vert%x(i)])
                ty = ConstructRealDynamicArray([streamlines(1)%y, vert%y(i)])

                ! Add
                call AddTopologicalMeshFace(topomesh, [thistp, i], &
                    tx, ty, TMfaceradID, 0_I8, 0.0_R8)
            end if 
        end do

        ! Housekeeping
        end associate

    end subroutine

#ifdef debug
    ! (X-point) streamlines
    subroutine AddTopologicalMeshStreamlines(topomesh, magneticField, &
        streamlinetracer, options)

        ! Description
        !============
        ! This routine adds any desired streamlines to the topological
        ! mesh. Currently, only streamlines originating from X-points 
        ! are implemented, as these are the only ones relevant for now.
        ! For X-points, we only trace the streamlines that do not go 
        ! to the core (these should already have been added in the
        ! extremum tracing routine). Furthermore, these streamlines are
        ! cut off at the first intersection with a poloidal face, or 
        ! they are not added if they first intersect with a non-poloidal
        ! face.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        type(magneticFieldUDT), intent(in)      :: magneticField
        class(StreamlineTracerUDT), intent(in)  :: streamlinetracer 
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: thisf
        integer(I8), allocatable, dimension(:)  :: tf, contourtypes, &
            fsIDs, tfnb
        real(R8)                                :: thisfval, traceval
        real(R8), allocatable, dimension(:)     :: xint, yint, tx, ty, &
            tfval, xtrace, ytrace
        logical                                 :: skipvertex
        logical, allocatable, dimension(:)      :: keepind, hasfv
        type(ContourUDT), allocatable           :: tc(:), allc(:)
        type(PolygonUDT)                        :: tcp 
        ! real(R8), allocatable, dimension(:)     :: 

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !============
        ! Associate
        associate(&
            nfs             => topomesh%nfs,    &
            face            => topomesh%face,   &
            vert            => topomesh%vert    &
        )

        ! Initialize
        allocate(allc(0), contourtypes(0), fsIDs(0))

        ! X-points
        !=========
        if (options%addxpstreamlines) then 
            ! Loop over all vertices, check for saddle points
            do i = 1, nfs 
                ! Check if the current vertex is a saddle point
                if (.not. (vert%type(i) == TMvertexsaddleID)) then 
                    cycle
                end if 

                ! For this saddle point, check which faces it has
                tf = vert%GetFace(i)
                
                ! Loop over all faces and find face neighbours. If the
                ! left neighbour and current face are both separatrix
                ! parts, trace
                do j = 1, size(tf) 
                    ! Check
                    if (.not. face%type(tf(j)) == TMfacesepID) then 
                        cycle
                    end if

                    ! Get neighbours
                    tfnb = vert%GetFaceNeig(i, tf(j))

                    ! Sanity check
                    if (size(tfnb) /= 2) then 
                        print *, 'vertex number: ', i
                        call gdErrorHandler('AddTopologicalMeshStreamlines: ' // & 
                            'could not find two face neighbours for vertex')
                    end if 

                    ! Check next face
                    if (.not. face%type(tfnb(1)) == TMfacesepID) then 
                        cycle
                    end if 

                    ! Compute starting point
                end do 

                ! Check how many 
            end do 
        end if 

        ! End association
        end associate

        ! Add to the topology mesh
        !=========================
        ! Add contours
        do  i = 1, size(allc)
            call InsertTopologicalMeshContour(topomesh, magneticField, &
                allc(i), contourtypes(i), fsIDs(i))
        end do 

        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Split boundaries
        call SplitTopologicalMeshFaces(topomesh)

    end subroutine
#endif 

    ! General region removal
    subroutine RemoveTopologicalMeshRegions(topomesh, vessel, options)

        ! Description
        !============
        ! Wrapper for all topological mesh region removal operations.
        ! Includes removal operations based on topology (e.g. core, 
        ! wide grid), vessel structure (based on vessel structure IDs), 
        ! and user input (simply indicated by user). See dedicated 
        ! subroutines for details on algorithms and implementation.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh 
        type(TopomeshOptionsUDT), intent(in)    :: options 
        type(VesselUDT), intent(in)             :: vessel

        ! Topology-based operations
        !==========================
        ! Core rggions
        if (options%removecoreregions) then 
            call RemoveTopologicalMeshCoreRegions(topomesh)
        end if 
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftercore')
        end if 

        ! Wide grid regions
        if (options%removewidegridregions) then 
            call RemoveTopologicalMeshWideGridRegions(topomesh)
        end if 

        ! Non-core regions
        if (options%removenoncoreregions) then 
            call RemoveTopologicalMeshNonCoreRegions(topomesh)
        end if

        ! Structure based operations
        !===========================
        if (options%removevesselregions) then 
            call RemoveTopologicalMeshVesselRegions(topomesh, vessel, options)
        end if 
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftervessel')
        end if 


    end subroutine

    ! Core region removal
    subroutine RemoveTopologicalMeshCoreRegions(topomesh)

        ! Description
        !============
        ! This routine removes extrema in the topological mesh (vertices type 1 or
        ! 3) and removes all their faces and regions. This routine is useful if one wants to
        ! generate grids that do not extend to the extrema in the field, but only a
        ! certain fraction (if used in combination with
        ! AddTopologicalMeshCoreBoundaries). 

        ! This routine should be used after the topological mesh has been fully
        ! constructed. 

        ! Notes 
        !======

        ! Declare variables
        !==================
        ! Arguments 
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv, tf
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i, j

        ! Initialize
        !===========
        ! Mark vertices for deletion
        delv = (topomesh%vert%type == TMvertexmaxID) .or. &
            (topomesh%vert%type == TMvertexminID)

        ! Mark faces for deletion
        allocate(delf(topomesh%face%ntot))
        delf = .false. 
        do i = 1, topomesh%face%ntot
            if (any(delv(topomesh%face%vert(i, :)))) then 
                delf(i) = .true.
            end if 
        end do

        ! Mark cells for deletion & convert aligned (non-deleted) faces
        ! to type TMfacecoreID
        allocate(delc(topomesh%cell%ntot))
        delc = .false.
        do i = 1, topomesh%cell%ntot
            tv = GetTMCellVert(topomesh%cell, i)
            if (any(delv(tv))) then 
                delc(i) = .true.
                tf = GetTMCellFace(topomesh%cell, i)
                do j = 1, size(tf)
                    if (topomesh%face%type(tf(j)) == TMfacepolID) then 
                        topomesh%face%type(tf(j)) = TMfacecoreID
                    end if 
                end do
            end if 
        end do 

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    ! Wide grid region removal
    subroutine RemoveTopologicalMeshWideGridRegions(topomesh)

        ! Description
        !============
        ! This routine removes any regions that are classified as 
        ! 'wide grid' regions. In practice, this means any region that 
        ! does not have a separatrix boundary segment (note that if 
        ! core boundaries are inserted, also the core region will be 
        ! removed, if not removed already)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i, j

        ! Initialize
        !===========
        ! Mark cells, vertices, faces for deletion
        allocate(delc(topomesh%cell%ntot), delv(topomesh%vert%ntot), &
            delf(topomesh%face%ntot))
        delc = .true. ! default true
        delv = .true.
        delf = .true.

        do i = 1, topomesh%cell%ntot 
            ! Get faces
            tcf = topomesh%cell%GetFace(i)

            ! Check if there is a separatrix face, then keep cell
            if (any(topomesh%face%type(tcf) == TMfacesepID)) then 
                delc(i)     = .false. 
                delf(tcf)   = .false.
                tcv         = topomesh%cell%GetVert(i)
                delv(tcv)   = .false.

                ! Set poloidal faces here to SOL faces
                do j = 1, size(tcf)
                    if (topomesh%face%type(tcf(j)) == TMfacepolID) then 
                        topomesh%face%type(tcf(j)) = TMfaceSOLID
                    end if 
                end do
            end if 
        end do

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    ! Remove all-but core regions (for Anthony)
    subroutine RemoveTopologicalMeshNonCoreRegions(topomesh)

        ! Description
        !============
        ! This routine removes all regions that don't have
        ! a maximum/minimum as boundary.
        ! This should yield the desired output for devices such as 
        ! e.g. TOMAS that are almost perfectly circular

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)              :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: delv, delf, delc

        ! Loop
        integer(I8)                         :: i, j

        ! Initialize
        !===========
        ! Mark cells, vertices, faces for deletion
        allocate(delc(topomesh%cell%ntot), delv(topomesh%vert%ntot), &
            delf(topomesh%face%ntot))
        delc = .true. ! default true
        delv = .true.
        delf = .true.

        do i = 1, topomesh%cell%ntot 
            ! Get faces
            tcf = topomesh%cell%GetFace(i)

            ! Get vertices
            tcv = topomesh%cell%GetVert(i)

            ! Check if there is an extremum, then keep cell and mark 
            ! aligned faces as core faces
            if (any(topomesh%vert%type(tcv) == TMvertexminID) .or. &
                any(topomesh%vert%type(tcv) == TMvertexmaxID)) then 
                delc(i)     = .false. 
                delf(tcf)   = .false.
                delv(tcv)   = .false.

                ! Check for aligned faces, set to core boundary
                do j = 1, size(tcf)
                    if (topomesh%face%type(tcf(j)) == TMfacepolID) then 
                        topomesh%face%type(tcf(j)) = TMfacecoreID 
                    end if 
                end do
            end if 
        end do

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    ! Remove regions that connect to some vessel structure(s)
    subroutine RemoveTopologicalMeshVesselRegions(topomesh, vessel, options)

        ! Description
        !============
        ! This routine removes (or retains, see later) any regions (
        ! topological mesh tubes) that has a radial boundary that 
        ! contains part of a vessel structure defined by the user. 
        ! The user can choose whether a radial line should be partly or 
        ! fully. The first one leads to most tubes deleted, the second 
        ! is more conservative. The user can also do the complementary 
        ! operations by instead keeping only these tubes that have this
        ! structure (partly) as a radial bound. Furthermore, if desired, 
        ! the user can request to cascade the deletion of tubes in the 
        ! radial direction either to increasing or decreasing psi values.
        ! This can be handy in cases with many different structures to 
        ! enumerate and where only a narrow grid is desired. 

        ! Algorithm
        !==========
        ! We use the vessel structure number interpolant to compute 
        ! the structure number(s) that each vessel topological face has.
        ! Then we check for each tube the radial faces. If they belong 
        ! (completely) to one or more specified structures, then they 
        ! are either marked for deletion or retainal. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh 
        type(VesselUDT), intent(in)             :: vessel 
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        logical                                 :: hasIDs
        logical, allocatable, dimension(:)      :: mark, &
            tempmark, newmark, delc, delf, delv
        integer(I8)                             :: vind
        integer(I8), allocatable, dimension(:)  :: tfvesselIDs, tf, &
            tc, tcf, tcv
        integer(I8), allocatable, dimension(:, :)   :: temp
        real(R8), allocatable, dimension(:)     :: xfv, yfv, xfe, yfe
        type(IntegerDynamicArrayUDT), allocatable   :: fvesselIDs(:)

        ! Loop
        integer(I8)                             :: i, j, k 

        ! Initialize
        !===========
        ! Unpack and associate
        associate(&
            tube        => topomesh%tube,   &
            face        => topomesh%face,   &
            vesselIDs   => options%rvrvesselIDs,    &
            docascade   => options%rvrdocascade,    &
            cascadedir  => options%rvrcascadedir,   &
            onlyfullycovered    => options%rvrfullycovered,     &
            doretain    => options%rvrretain        &
            )

        ! Hedge for edge cases
        if (size(vesselIDs) == 0) then 
            ! Nothing to do here, move along
            return 
        end if 
        
        ! Mark vector
        allocate(mark(tube%ntot), fvesselIDs(face%ntot))
        mark = .false. 

        ! Determine face structure IDs
        !=============================
        do i = 1, face%ntot 
            ! Check if it is a structure
            if (face%type(i) == TMfacebndID) then 
                ! Get coordinates of face edge centers
                xfv = face%x(i)%Get()
                yfv = face%y(i)%Get()
                xfe = 0.5*(xfv(1:size(xfv)-1) + xfv(2:))
                yfe = 0.5*(yfv(1:size(yfv)-1) + yfv(2:))

                ! Compute vessel IDs
                call vessel%exactplfvessel%EvaluateLabel(xfe, yfe, temp)
                call Unique(temp(:, 1), tfvesselIDs)

                ! Assign
                fvesselIDs(i) = ConstructIntegerDynamicArray(tfvesselIDs)
            else
                fvesselIDs(i) = ConstructIntegerDynamicArray()
            end if 
        end do 

        ! Mark tubes
        !===========
        ! Mark tubes based on structure IDs
        do i = 1, tube%ntot
            ! Get tube radial faces
            tf = GetTMTubeFace(tube, i)

            ! Compare IDs with user-defined IDs
            do j = 1, size(tf)
                ! Get structure IDs
                tfvesselIDs = fvesselIDs(tf(j))%Get()

                ! Skip if no vessel IDs are present
                if (size(tfvesselIDs) == 0) then 
                    cycle
                end if 

                ! Compare
                if (onlyfullycovered) then 
                    hasIDs = .true. 
                    do k = 1, size(tfvesselIDs)
                        hasIDs = hasIDs .and. any(tfvesselIDs(k) == vesselIDs)
                    end do 
                    if (hasIDs) then 
                        mark(i) = .true. 
                    end if 
                else
                    hasIDs = .false.
                    do k = 1, size(tfvesselIDs)
                        hasIDs = hasIDs .or. any(tfvesselIDs(k) == vesselIDs)
                    end do 
                    if (hasIDs) then 
                        mark(i) = .true. 
                    end if 
                end if 

                ! If marked -> exit
                if (mark(i)) then 
                    exit 
                end if
            end do 
        end do 

        ! Mark additional tubes by cascade if desired
        if (docascade) then 

            ! Initialize
            newmark = mark 
            newmark = .false. 
            
            select case (cascadedir)

            case ('none', 'no')

                ! Do nothing

            case ('upstream', 'up')

                ! Need to flood from high to low -> need to reverse direction
                ! of graph
                call tube%graph%ReverseDirection()

            case ('downstream', 'down')

                ! Need to flood from low to high -> can just flood as is


            case default

                call gdErrorHandler('RemoveTopologicalMeshVesselRegions: ' // & 
                    'unknown cascade direction: ' // cascadedir)

            end select

            ! Mark tubes
            do i = 1, tube%ntot 
                if (mark(i)) then 
                    vind = tube%graph%GetVertexIndex(i)
                    tempmark = tube%graph%Flood(vind)
                    newmark = newmark .or. tempmark
                end if 
            end do 
            mark = mark .or. newmark

            ! Check if we need to re-reverse
            select case (cascadedir)

            case ('upstream')

                ! Re-reverse to get original direction again
                call tube%graph%ReverseDirection()

            end select ! exception handling done before

        end if 

        ! Check if we need to retain instead of delete
        if (doretain) then 
            mark = .not. mark 
        end if 

        ! Housekeeping
        end associate

        ! Delete
        !=======
        ! Initialize
        !-----------
        ! Mark cells, vertices, faces for deletion - note: we actually
        ! work the other way around, i.e. we only keep tubes that are 
        ! not to be deleted. This yields immediately also the correct 
        ! faces, cells, and vertices to keep. 
        allocate(delc(topomesh%cell%ntot), delv(topomesh%vert%ntot), &
            delf(topomesh%face%ntot))
        delc = .true. ! default true
        delv = .true.
        delf = .true.

        do i = 1, topomesh%tube%ntot 
            ! Get cells, faces, vertices and mark for retainal
            tc = GetTMTubeCell(topomesh%tube, i)

            ! Check if tube should be deleted
            if (.not. mark(i)) then ! This tube should be retained
                delc(tc) = .false. 
                do j = 1, size(tc)
                    tcf = GetTMCellFace(topomesh%cell, tc(j))
                    tcv = GetTMCellVert(topomesh%cell, tc(j))
                    delf(tcf) = .false. 
                    delv(tcv) = .false. 
                end do 
            else
                ! This tube should be deleted - set aligned poloidal 
                ! face types to SOL (note: this may be wrong in some
                ! exceptional cases if tubes are deleted by cascading adn
                ! haphazardly some core regions still remain)
                delc(tc) = .true. 
                do j = 1, size(tc)
                    tcf = GetTMCellFace(topomesh%cell, tc(j))
                    do k = 1, size(tcf)
                        if (topomesh%face%type(tcf(k)) == TMfacepolID) then 
                            topomesh%face%type(tcf(k)) = TMfaceSOLID
                        end if 
                    end do 
                end do 
            end if 
        end do

        ! Delete
        !=======
        ! Vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, delv)

        ! Faces
        call RemoveTopologicalMeshFaceLogical(topomesh, delf)

        ! Cells
        call RemoveTopologicalMeshCellLogical(topomesh, delc)

        ! Update again
        !=============
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    !------------------------------------------------------------------!
    !                       TOPOLOGICAL MESH ADAPTOR                   !
    !------------------------------------------------------------------!

    ! Topomesh adaptor initialization
    subroutine InitializeTopomeshAdaptor(tmadaptor, topomesh, &
        fieldtracer, streamlinetracer, magneticField, vessel, options)

        ! Description
        !============
        ! This subroutine initializes the topological mesh adaptor

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor 
        type(TopomeshUDT), intent(in)           :: topomesh 
        class(ContourtracerUDT), intent(in)     :: fieldtracer
        class(StreamlinetracerUDT), intent(in)  :: streamlinetracer
        type(magneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(in)             :: vessel
        type(TopomeshOptionsUDT), intent(in)    :: options 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: temp, temp_orig

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Associate
        associate(&
            face    => topomesh%face,   &
            tube    => topomesh%tube    &
            )

        ! Check allocation status
        if (allocated(tmadaptor%illegalfsIDs)) deallocate(tmadaptor%illegalfsIDs)
        if (allocated(tmadaptor%splitfsIDs)) deallocate(tmadaptor%splitfsIDs)
        if (allocated(tmadaptor%facepsi)) deallocate(tmadaptor%facepsi)
        if (allocated(tmadaptor%facepsi_mono)) deallocate(tmadaptor%facepsi_mono)
        if (allocated(tmadaptor%facedlcrad)) deallocate(tmadaptor%facedlcrad)
        if (allocated(tmadaptor%tubedpsi)) deallocate(tmadaptor%tubedpsi)
        if (allocated(tmadaptor%tubelrad)) deallocate(tmadaptor%tubelrad)
        if (allocated(tmadaptor%fieldtracer)) deallocate(tmadaptor%fieldtracer)
        if (allocated(tmadaptor%streamlinetracer)) deallocate(tmadaptor%streamlinetracer)

        ! Initialize
        allocate(tmadaptor%illegalfsIDs(0), tmadaptor%facepsi(face%ntot), &
            tmadaptor%facedlcrad(face%ntot), tmadaptor%tubedpsi(tube%ntot), &
            tmadaptor%tubelrad(tube%ntot), tmadaptor%splitfsIDs(0), &
            tmadaptor%facepsi_mono(face%ntot))
        allocate(tmadaptor%fieldtracer, source=fieldtracer)

        ! Copy (for face data updating later on)
        tmadaptor%fieldtracer   = fieldtracer
        tmadaptor%streamlinetracer  = streamlinetracer
        tmadaptor%magneticField = magneticField
        tmadaptor%vessel        = vessel

        ! Copy options
        tmadaptor%dpsimin = options%dpsimintubes
        tmadaptor%lradmin = options%lradmintubes
        tmadaptor%lpolsepmin = options%lpolsepmin
        tmadaptor%allowsepmerge = options%mtallowseparatrix
        tmadaptor%allowsepsimplify = options%simplifyseparatrix
        tmadaptor%allowcoremerge = options%mtallowcore
        tmadaptor%allowpfmerge  = options%mtallowpf
        tmadaptor%mergetubemeth = options%mergetubemeth

        ! Compute data
        !=============
        call wall_time(ts)
        ! Face psi values and radial length
        !$omp parallel do default(none) if (.not. omp_in_parallel()) &
        !$omp private(i, temp, temp_orig) &
        !$omp shared(topomesh, fieldtracer, tmadaptor, magneticField)
        do i = 1, face%ntot
            ! Compute
            call GetTMFacePsiValueDistribution(topomesh, fieldtracer, i, &
                temp, temp_orig)
            tmadaptor%facepsi(i) = ConstructRealDynamicArray(temp_orig)
            tmadaptor%facepsi_mono(i) = ConstructRealDynamicArray(temp)
            temp = GetTMFaceRadialLengthDistribution(topomesh, &
                fieldtracer, magneticField, i)
            tmadaptor%facedlcrad(i) = ConstructRealDynamicArray(temp)
        end do 
        !$omp end parallel do
        call wall_time(te)
        print *, 'time spent in adaptor setup = ', te - ts

        ! Housekeeping
        end associate

    end subroutine

    ! Merge criterion evaluation
    subroutine EvaluateTMTubesMergeCriterionTA(tmadaptor, topomesh, &
        includealbndin)

        ! Description
        !============
        ! This routine evaluates all possible merging criteria - currently
        ! only the minimal psi width and minimal radial length criteria.
        ! These values are computed for all tubes. 
        
        ! Note: we now also allow passing the 'includealbndin' logical
        ! which will lead to inclusion of aligned boundary faces if 
        ! set to true. IF not present, default value is false. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(in)           :: topomesh 
        logical, optional, intent(in)           :: includealbndin

        ! Auxiliary
        real(R8)                                :: lrad, psimin, &
            psimax
        integer(I8), allocatable, dimension(:)  :: tf
        logical                                 :: includealbnd

        ! Loop
        integer(I8)                     :: i

        ! Initialize
        !===========
        ! Check input
        if (present(includealbndin)) then 
            includealbnd = includealbndin
        else
            includealbnd = .false.
        end if

        ! Check sizes
        if ((size(tmadaptor%tubedpsi) /= topomesh%tube%ntot) .or. &
            (size(tmadaptor%tubelrad) /= topomesh%tube%ntot))then 
            deallocate(tmadaptor%tubedpsi, tmadaptor%tubelrad)
            allocate(tmadaptor%tubedpsi(topomesh%tube%ntot), &
                tmadaptor%tubelrad(topomesh%tube%ntot))
        end if 

        ! Compute
        !========
        do i = 1, topomesh%tube%ntot
            ! Get tube faces
            tf = topomesh%tube%GetFace(i)

            ! Evaluate psi criterion
            call GetTMTubePsiLimitsTA(tmadaptor, topomesh, i, psimin, psimax, includealbndin=includealbnd)
            tmadaptor%tubedpsi(i) = max(psimax - psimin, 0.0_R8) 

            ! Evaluate radial length criterion
            call GetTMTubeRadialWidthTA(tmadaptor, topomesh, &
                i, lrad, includealbndin=includealbnd)
            tmadaptor%tubelrad(i) = lrad

        end do 

    end subroutine

    ! Flux tube merge driver
    subroutine MergeTMTubesDriverTA(tmadaptor, topomesh, options)

        ! Description
        !============
        ! This routine is the main driver for tube merging. The main 
        ! idea is that first simple merges are performed (if possible) 
        ! before more complex merges involving multiple tubes are done.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        logical                             :: wasmerged, &
            appliedsplitting

        ! Initialize
        !===========
        ! Check if we should merge, otherwise return
        if (.not. options%mergetubes) then 
            return 
        end if
        wasmerged = .false.
        appliedsplitting = .false.

        ! Loop until all are merged
        do while (.true.) 
            ! Write output
            if (options%writedebugoutput) then 
                call WriteTopologicalMesh(topomesh, 'topomesh_temp')
            end if 

            ! Is there a simple merge that can be done? 
            if (.not. appliedsplitting) then ! skip if we splitted previously, already doing complex merge
                call tmadaptor%MergeTMTubesSimple(topomesh, wasmerged)
                if (wasmerged) then 
                    if (options%writedebugoutput) then 
                        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
                    end if 
                    print *, 'MergeTopologicalMeshFluxTubes: applied simple merge'
                    cycle
                end if 
            end if 

            ! Is there a complex merge that can be done?
            call tmadaptor%MergeTMTubesComplex(topomesh, options, wasmerged, &
                appliedsplitting)

            ! In case we merged a set of tubes, we can reset the split
            ! flux surface IDs (those may be merged over now)
            if (wasmerged .and. .not. appliedsplitting) then 
                deallocate(tmadaptor%splitfsIDs)
                allocate(tmadaptor%splitfsIDs(0))
            end if 

            ! Write output
            if (wasmerged .or. appliedsplitting) then 
                if (wasmerged) then 
                    if (options%writedebugoutput) then 
                        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
                    end if 
                    print *, 'MergeTopologicalMeshFluxTubes: applied complex merge'
                elseif (appliedsplitting) then 
                    if (options%writedebugoutput) then 
                        call WriteTopologicalMesh(topomesh, 'topomesh_temp')
                    end if 
                    print *, 'MergeTopologicalMeshFluxTubes: applied splitting'
                end if
                cycle
            end if 

            ! If no merging operations have been done, exit
            if (.not. (wasmerged .or. appliedsplitting)) then 
                exit 
            end if 
        end do 

        ! Clean up splitted tubes
        deallocate(tmadaptor%splitfsIDs)
        allocate(tmadaptor%splitfsIDs(0))
        call tmadaptor%CleanSplitTMTubes(topomesh)

    end subroutine

    ! Simple flux tube merging
    subroutine MergeTMTubesSimpleTA(tmadaptor, topomesh, wasmerged)

        ! Description
        !============
        ! This routine merges a single flux tube by removing its poloidal
        ! boundaries at one side, if possible. Only 'simple' merges are 
        ! allowed, meaning that the flux tube only can have one 
        ! neighbour, and that the radial faces at both sides of both 
        ! tubes connect with each other. This is the most basic merging
        ! operation, which will likely not occur in the topomesh since 
        ! these kind of tube setups are often avoided from the start. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh 
        logical, intent(out)                    :: wasmerged
        
        ! Auxiliary
        integer(I8)                             :: tubepairind, t1, t2 
        integer(I8), allocatable, dimension(:)  :: tf1, tf2, &
            tube1, tube2
        logical, allocatable, dimension(:)      :: ismarked, ismarkedpair
        real(R8), allocatable, dimension(:)     :: dpsi, dlrad, dpsipair, &
            dlradpair
        type(IntegerDynamicArrayUDT), allocatable   :: hftubes(:), lftubes(:)

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !==========
        ! Unpack
        associate(&
            face    => topomesh%face,   &
            tube    => topomesh%tube    &
            )

        ! Set output
        wasmerged = .false. 

        ! Compute criteria
        call tmadaptor%EvaluateTMTubesMergeCriterion(topomesh)
        dpsi = tmadaptor%tubedpsi - tmadaptor%dpsimin 
        dlrad = tmadaptor%tubelrad - tmadaptor%lradmin
        ismarked = (dpsi < 0.0_R8) .or. (dlrad < 0.0_R8)

        ! Easy preliminary check
        if (.not. any(ismarked)) then 
            return 
        end if 

        ! Determine mergeable tube pairs
        !===============================
        ! Call dedicated subroutine
        call tmadaptor%GetMergeTubePairs(topomesh, hftubes, lftubes, ismarked)

        ! Check which pairs are simple and are allowed to be merged
        allocate(tube1(0), tube2(0))
        do i = 1, size(hftubes)
            ! Check sizes
            if ((hftubes(i)%Size() /= 1) .or. (lftubes(i)%Size() /= 1)) then 
                cycle 
            end if 

            ! Unpack
            t1 = hftubes(i)%Get(1)
            t2 = lftubes(i)%Get(1)

            ! Check if the current tube is closed -> only allowed to merge
            ! with other closed tubes
            if (tube%isclosed(t1)) then 
                ! Check tube neighbours at first side
                if (tube%isclosed(t2)) then 
                    ! Tubes should be mergeable unless something is 
                    ! wrong with the initial topomesh
                    tube1 = [tube1, t1]
                    tube2 = [tube2, t2]
                end if 
            else
                ! Check tube neighbours at first side
                if (.not. tube%isclosed(t2)) then 
                    ! Check if start and end radial faces have vertices
                    ! in common
                    tf1 = tube%GetFace(t1)
                    tf2 = tube%GetFace(t2)
                    if ((AreTMFacesAdjacent(face, tf1(1), tf2(1)) .and. &
                        AreTMFacesAdjacent(face, tf1(size(tf1)), tf2(size(tf2)))) .or. &
                        (AreTMFacesAdjacent(face, tf1(1), tf2(size(tf2))) .and. &
                        AreTMFacesAdjacent(face, tf1(size(tf1)), tf2(1)))) then 
                        
                        ! Tubes should be mergeable
                        tube1 = [tube1, t1]
                        tube2 = [tube2, t2]
                    end if 
                end if 
            end if 
        end do 

        ! Housekeeping
        end associate

        ! Merge tubes
        !============
        ! Check which pairs are marked (should be all actually)
        ismarkedpair = ismarked(tube1) .and. ismarked(tube2)

        ! Set illegal pairs to be unmarked
        do i = 1, size(ismarkedpair)
            ismarkedpair(i) = ismarkedpair(i) .and. &
                tmadaptor%IsTubePairMergeable(topomesh, [tube1(i)], [tube2(i)])
        end do 

        ! Give priority to 'smallest' tube pair (in terms of psi difference)
        dpsipair = dpsi(tube1) + dpsi(tube2)
        dlradpair = dlrad(tube1) + dlrad(tube2)
        where (.not. ismarkedpair) 
            dpsipair = posinfval_R8()
            dlradpair = posinfval_R8()
        end where
        
        ! Loop until merged or no tube pairs available for merge anymore
        do while (any(ismarkedpair) .and. .not. wasmerged)
            ! Find the next possible pair
            tubepairind = minloc(dpsipair, 1)

            ! Set values to posinf to ignore in next loop if necessary
            dpsipair(tubepairind) = posinfval_R8()
            dlradpair(tubepairind) = posinfval_R8()
            ismarkedpair(tubepairind) = .false. 

            ! Attempt to merge
            call tmadaptor%MergeTMTubes(topomesh, [tube1(tubepairind)], [tube2(tubepairind)], &
                wasmerged)

        end do 

    end subroutine

    ! Complex flux tube pair merging
    subroutine MergeTMTubesComplexTA(tmadaptor, topomesh, options, &
        wasmerged, appliedsplitting)

        ! Description
        !============
        ! This routine performs more complicated merges between multiple
        ! flux tubes that may be open or closed. To minimize 
        ! impact on 'proper' flux tubes, additional contours may be
        ! inserted to first split existing tubes into one small and one
        ! large tube, after which the smaller tube is merged with the 
        ! original small tube. 

        ! Algorithm
        !==========
        ! 1) Determine which tubes are 'too small' based on merge criterion
        ! 2) Determine tube sets that may be merged for these tubes. If 
        ! one tube set is found that contains only too small tubes, merge 
        ! this tubeset. If multiple tube sets are found, then prioritize
        ! merging of the 'smallest' tube set (in terms of delta psi for example)
        ! Otherwise, go to 3)
        ! 3) In this case, all remaining tube sets contain at least one
        ! tube that is considered too large to merge. Check for each tube
        ! if we can insert a contour to split this tube into two (or three, 
        ! if on the other side there's also a tube to be merged) tubes.
        ! The smaller tube should be that small such that it would be 
        ! marked for merging. 
        ! 4) Insert contours where applicable into the topomesh and 
        ! rebuild. Go to 1)

        ! Notes
        !======
        ! Note 1: in practice, we are agnostic that this routine is called
        ! by an overarching routine that loops until all desired/possible
        ! merges are done. Therefore, in the case of step 3, we only 
        ! apply the splitting operation and then return, setting 
        ! 'appliedsplitting' to true (and wasmerged to false). In this 
        ! case, the overarching routine can recall this function to 
        ! perform the merge, and we avoid having to make this a 
        ! recursive subroutine with possibly undesired side effects

        ! Note 2: we prevent for now merging of tubes that would result
        ! in a final tube with overlapping psi values. In many 
        ! conventional setups, this would likely not happen, but for more
        ! complex vessel geometries this is actually a potentially common 
        ! case. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh 
        type(TopomeshOptionsUDT), intent(in)    :: options 
        logical, intent(out)                    :: wasmerged, &
            appliedsplitting

        ! Auxiliary
        integer(I8)                             :: tubepairind
        integer(I8), allocatable, dimension(:)  :: tube1, tube2, &
            splittubes
        logical, allocatable, dimension(:)      :: ismarked, &
            ismarkedpair, issplittable, dolfside, dohfside, ismergeable
        real(R8), allocatable, dimension(:)     :: dvalpair, dpsi, dlrad
        type(IntegerDynamicArrayUDT)            :: splittubesida
        type(IntegerDynamicArrayUDT), allocatable   :: hftubes(:), &
            lftubes(:)

        ! Loop
        integer(I8)                             :: i, tpc

        ! Initialize
        !===========
        ! Set output
        wasmerged = .false. 
        appliedsplitting = .false.

        ! Compute criteria 
        ! Note: no aligned boundaries are included here because tubes 
        ! may have overlapping psi values because of that, at which point they're not 
        ! mergeable anyway anymore...
        call tmadaptor%EvaluateTMTubesMergeCriterion(topomesh, includealbndin=.false.) ! don't include aligned boundary faces here 
        dpsi = tmadaptor%tubedpsi - tmadaptor%dpsimin 
        dlrad = tmadaptor%tubelrad - tmadaptor%lradmin
        ismarked = (dpsi < 0.0_R8) .or. (dlrad < 0.0_R8)

        ! Easy preliminary check
        if (.not. any(ismarked)) then 
            return 
        end if 

        ! Check which tubes are splittable (if not splittable, it is 
        ! mergeable but not necessarily marked)
        issplittable = tmadaptor%IsTubeSplittable(topomesh)

        ! Determine mergeable tube pairs
        !===============================
        ! Call dedicated subroutine
        call tmadaptor%GetMergeTubePairs(topomesh, hftubes, lftubes, ismarked)
        tpc = size(hftubes)

        ! If no pairs were found, exit
        if (tpc == 0) then 
            print *, 'MergeTopologicalMeshFluxTubesComplex: no tube pairs ' // &
                'found that could be merged, exiting...'
            return 
        end if 

        ! Merge tubes
        !============
        ! Check if there are any tube pairs that consist fully of 
        ! mergeable pairs (i.e. none are splittable)
        allocate(ismarkedpair(tpc), dvalpair(tpc), ismergeable(tpc))
        ismarkedpair = .false. 
        ismergeable = .false. 
        do i = 1, tpc
            ! Get tubes
            tube1 = hftubes(i)%Get()
            tube2 = lftubes(i)%Get()

            ! Check if tube pair is legal
            ismergeable(i) = tmadaptor%IsTubePairMergeable(topomesh, tube1, tube2)
            if (.not. ismergeable(i)) then 
                ismarkedpair(i) = .false.
                dvalpair(i) = posinfval_R8()
                cycle
            end if 

            ! Check if all tubes are mergeable if they are legal
            if (.not. any(issplittable([tube1, tube2]))) then 
                ismarkedpair(i) = .true.
                dvalpair(i) = minval(dpsi(tube1)) + minval(dpsi(tube2))
            else
                ismarkedpair(i) = .false.
                dvalpair(i) = posinfval_R8()
            end if 
        end do

        ! If any pairs were found, merge these
        do while (any(ismarkedpair) .and. .not. wasmerged)
            ! Find the next possible pair
            tubepairind = minloc(dvalpair, 1)

            ! Set values to posinf to ignore in next loop if necessary
            dvalpair(tubepairind) = posinfval_R8()
            ismarkedpair(tubepairind) = .false. 

            ! Attempt to merge
            call tmadaptor%MergeTMTubes(topomesh, hftubes(tubepairind)%Get(), lftubes(tubepairind)%Get(), &
                wasmerged)

        end do 

        ! Split tubes
        !============
        ! If no tubes were merged, check if we can/should split tubes
        if (.not. wasmerged) then 
            ! Initialize
            splittubesida = ConstructIntegerDynamicArray()
            allocate(dolfside(topomesh%tube%ntot))
            dolfside = .false. 
            dohfside = dolfside

            ! Check if we should apply any splitting
            do i = 1, tpc
                if (ismergeable(i)) then ! Don't allow splitting for non-mergeable - will lead to cyclic behavior of merge/split!
                    ! Get tubes
                    tube1 = hftubes(i)%Get()
                    tube2 = lftubes(i)%Get()
                    
                    ! Check which side to split
                    if (any(ismarked(tube1)) .or. any(ismarked(tube2))) then ! tube 1 is hf side
                        dohfside(tube2) = .true. ! Need to split at high flux side of tube 2
                        dolfside(tube1) = .true. ! Need to split at low field side of tube 1 
                        call splittubesida%Append(pack(tube2, issplittable(tube2)))
                        call splittubesida%Append(pack(tube1, issplittable(tube1)))
                    end if 
                end if 
            end do

            ! Split tubes based on merge criterion and side - make sure only unique tubes are retained
            call Unique(splittubesida%Get(), splittubes)
            if (size(splittubes) > 0) then 
                appliedsplitting = .true.
                call tmadaptor%SplitTMTubes(topomesh, splittubes, dohfside(splittubes), &
                    dolfside(splittubes), options)
            end if 
        end if 

    end subroutine

    ! Topological tube merging operator 
    subroutine MergeTMTubesTA(tmadaptor, topomesh, hftubes, lftubes, &
        wasmerged)

        ! Description
        !============
        ! This routine merges two sets of tubes given by their IDs in
        ! hftubes and lftubes by removing the exising aligned faces 
        ! between them and introducing aligned vessel parts where 
        ! needed. It is assumed that it is checked beforehand whether 
        ! these tubes are allowed to merge or not. Additionally, it is 
        ! assumed that the common faces form a single (possibly branching)
        ! polygon. If multiple are present, an error is currently thrown
        ! (we could also just ignore the other tubes but OK). Given that
        ! we can rather easily know beforehand which tubes are the high
        ! flux and low flux side, we assume this is given before. 

        ! Algorithms
        !===========
        ! Main algorithm
        !---------------
        ! 1) Find common aligned tube faces to determine which boundary
        ! will be merged. 
        ! 2) Sort these aligned faces 
        ! 3) Check if the sorted faces form a closed polygon or not. If
        ! they do not, then all tubes should be open. Otherwise, at least
        ! one tube should be closed (and should be the only tube on one
        ! side normally speaking). If it's a branching polygon, then it
        ! should be a separatrix. 
        ! 4) Depending on 3), apply either:
        !   - open-open merge (both sides have only open tubes)
        !   - open-closed merge (one side has a closed tube but the merging faces do not belong to a separatrix)
        !   - closed-closed merge (both sides have closed tubes)
        !   - separatrix merge (face polygon is branching) (not supported currently)
        

        ! open-open merge
        !----------------
        ! 1) Find the radial faces at the boundary aligned faces (i.e. 
        ! the first and the last one) of both tubes. These faces should
        ! be adjacent to the boundary vertex of the boundary aligned 
        ! face (at least one). Find all boundary faces that are adjacent
        ! to these without having an aligned face inbetween - all these
        ! boundary faces kan be kept as actual vessel faces.
        ! 2) Retype all other vessel boundary faces that occur in the 
        ! tubes to be aligned vessel parts (this in order to achieve 
        ! again a continuous flux surface boundary at both sides of the
        ! new tube that will originate from the merge)
        ! 3) Remove all non-boundary topological mesh flux surfaces that
        ! were marked before
        ! 4) Simplify the topomesh to remove any remaining unecessary 
        ! points and reconstruct all topomesh quantities

        ! closed-closed merge
        !--------------------
        ! This can only happen for non-branching polygons for flux 
        ! surfaces that are e.g. core surfaces. 
        ! 1) Determine which radial faces of both tubes are adjacent.
        ! These faces can be kept in the topomesh, others will be 
        ! deleted
        ! 2) Remove the previously marked faces and the common aligned
        ! faces. Normally, in a closed-closed merge, there cannot be
        ! any boundary faces or aligned vessel parts
        ! 3) Simplify the topomesh to remove any remaining unecessary 
        ! points and reconstruct all topomesh quantities

        ! open-closed merge
        !------------------
        ! This should only happen for limiter cases if the faces do 
        ! not form a branching polygon. In that case, the radial faces
        ! of the open tube should be adjacent to the radial face of the
        ! closed tube. Merging is done by removing the limiter curve and
        ! retyping the radial faces of the open tube as aligned vessel
        ! parts. 

        ! separatrix merge
        !-----------------
        ! This is currently not implemented/supported for several 
        ! reasons:
        ! - it is very tricky to do since regions will inevitably be 
        ! destroyed
        ! - it is very often not desired: removing these regions will
        ! lead to very ugly grids and probably very difficult to converge
        ! simulations
        ! - the only use case is probably to remove small scale behavior
        ! of the magnetic field. In that case, it is probably better
        ! to reduce the magnetic field resolution in order to filter 
        ! these effects out, or to use the existing adjustments on the 
        ! contour level to e.g. allow connected null configurations. 

        ! Should this ever be implemented, one can consider the following
        ! merge cases:
        ! - full deletion of an X-point by recursively deleting all 
        ! connecting regions upstream or downstream to the separatrix 
        ! part
        ! - merging of nearly coinciding X-points into a connected 
        ! configuration 
        ! - deletion of PF-like region close to the wall
        ! This is obviously not an exhaustive list, but should cover 
        ! some rather plausible use cases. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)                   :: tmadaptor
        type(TopomeshUDT), intent(inout)            :: topomesh
        integer(I8), dimension(:), intent(in)       :: hftubes, lftubes
        logical, intent(out)                        :: wasmerged

        ! Auxiliary
        integer(I8)                                 :: ne
        integer(I8), allocatable, dimension(:)      :: hfface1, &
            hfface2, lfface1, lfface2, tf, tf1, tf2, sortind, polygonID, &
            mergefaces, mergevert
        integer(I8), allocatable, dimension(:, :)   :: tfv
        logical                                     :: isbranchingmerge
        logical, allocatable, dimension(:)          :: ispolygonstart, &
            isbranchingpolygon

        ! Loop
        integer(I8)                                 :: i

        ! Initialize
        !===========
        ! Initialize
        wasmerged = .false. 

        ! Unpack
        associate(&
            tube    => topomesh%tube,   &
            face    => topomesh%face    &
            )

        ! Sanity checks
        if ((size(hftubes) == 0) .or. (size(lftubes) == 0)) then 
            ! This shouldn't happen, simply return
            print *, 'MergeTMTubes: no tubes at at least one side, ' // & 
                'no merge to be done. Returning...'
            return 
        end if 

        ! Get all high and low field faces of tubes
        allocate(hfface1(0), hfface2(0), lfface1(0), lfface2(0))
        do i = 1, size(hftubes)
            tf1 = tube%GetHighFluxBndFace(hftubes(i))
            tf2 = tube%GetLowFluxBndFace(hftubes(i))
            hfface1 = [hfface1, tf1]
            lfface1 = [lfface1, tf2]
        end do 
        do i = 1, size(lftubes)
            tf1 = tube%GetHighFluxBndFace(lftubes(i))
            tf2 = tube%GetLowFluxBndFace(lftubes(i))
            hfface2 = [hfface2, tf1]
            lfface2 = [lfface2, tf2]
        end do 

        ! Sanity checks
        tf = GetCommonElements(hfface1, hfface2)
        if (size(tf) > 0) then 
            ! Both tubes have high field faces in common, this is unexpected
            call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
            call gdErrorHandler('MergeTMTubes: tubes have high field ' // & 
                'faces in common, unexpected')
        end if 
        tf = GetCommonElements(lfface1, lfface2)
        if (size(tf) > 0) then 
            ! Both tubes have low field faces in common, this is unexpected
            call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
            call gdErrorHandler('MergeTMTubes: tubes have low field ' // & 
                'faces in common, unexpected')
        end if 

        ! Get common faces
        tf1 = GetCommonElements(hfface1, lfface2)
        tf2 = GetCommonElements(hfface2, lfface1)

        ! Sanity checks
        if ((size(tf1) > 0) .and. (size(tf2) > 0)) then 
            ! Both tubes have low field faces in common, this is unexpected
            call gdErrorHandler('MergeTMTubes: tubes have common ' // & 
                'faces at both sides, unexpected')
        elseif ((size(tf1) == 0) .and. (size(tf2) == 0)) then 
            ! No common faces, also not expected
            call gdErrorHandler('MergeTMTubes: tubes have no faces in common, ' // & 
                'unexpected')
        elseif (size(tf1) > 0) then 
            ! Make sure to take all faces and not only common faces! Otherwise, 
            ! some boundary faces that are crucial for the merging operation
            ! will not be considered...
            call Unique([hfface1, lfface2], tf)
            !tf = GetCommonElements(hfface1, lfface2)
        else 
            call Unique([hfface2, lfface1], tf)
            !tf = GetCommonElements(hfface2, lfface1)
        end if 

        ! Sort faces by vertices
        tfv = face%vert(tf, :)
        ne = size(tf)
        allocate(sortind(ne), ispolygonstart(ne), isbranchingpolygon(ne))
        call SortPolygonEdges(tfv, ne, sortind, &
            ispolygonstart, isbranchingpolygon, polygonID)
        tfv(:, 1) = tfv(sortind, 1)
        tfv(:, 2) = tfv(sortind, 2)
        tf = tf(sortind)

        ! Check if only one polygon (possibly branching) is found
        if (maxval(polygonID) > 1) then 
            call gdErrorHandler('MergeTMTubes: multiple polygons found, ' // & 
                'probably too many tubes passed for merging. ')
        end if 

        ! Extract sorted polygon faces (note: for branching polygons, 
        ! all faces are simply appended now)
        mergefaces = tf 
        isbranchingmerge = .false. 
        if (any(isbranchingpolygon)) then 
            isbranchingmerge = .true.
        end if

        ! Housekeeping
        end associate

        ! Check the merge case
        !=====================
        if (isbranchingmerge) then 
            ! Separatrix merge 
            call tmadaptor%MergeTMTubesS(topomesh, hftubes, lftubes, &
                mergefaces, wasmerged)

            ! Issue message in case of failure
            if (.not. wasmerged) then 
                print *, 'MergeTMTubesTA: separatrix merge not successful, continuing...'
            end if 
        else
            ! Check faces and vertices
            call ExtractPolygonVertices(topomesh%face%vert(mergefaces, :), &
                size(mergefaces), mergevert)

            ! Check for closed face
            if (mergevert(1) == mergevert(size(mergevert))) then 
                ! Closed polygon, check tubes

                ! Some sanity checks, only two tubes expected
                if ((size(hftubes) /= 1) .or. size(lftubes) /= 1) then 
                    call gdErrorHandler('MergeTMTubesTA: expected ' // & 
                        'only one tube on each side for closed-closed merge, ' // & 
                        'but obtained multiple ones - unexpected.')
                end if 

                ! Check merge type
                if (any(topomesh%tube%isclosed(hftubes)) .and. &
                    any(topomesh%tube%isclosed(lftubes))) then 

                    ! Closed-closed merge
                    call tmadaptor%MergeTMTubesCC(topomesh, hftubes(1), lftubes(1), &
                        mergefaces, wasmerged)

                    ! Issue message in case of failure
                    if (.not. wasmerged) then 
                        print *, 'MergeTMTubesTA: closed-closed merge not successful, continuing...'
                    end if 

                elseif (any(topomesh%tube%isclosed(hftubes))) then 

                    ! Open-closed merge
                    call tmadaptor%MergeTMTubesOC(topomesh, hftubes(1), lftubes(1), &
                        mergefaces, wasmerged)

                    ! Issue message in case of failure
                    if (.not. wasmerged) then 
                        print *, 'MergeTMTubesTA: open-closed merge not successful, continuing...'
                    end if 

                elseif (any(topomesh%tube%isclosed(lftubes))) then 

                    ! Open-closed merge
                    call tmadaptor%MergeTMTubesOC(topomesh, hftubes(1), lftubes(1), &
                        mergefaces, wasmerged)

                    ! Issue message in case of failure
                    if (.not. wasmerged) then 
                        print *, 'MergeTMTubes: open-closed merge not successful, continuing...'
                    end if 

                else 

                    ! Open-open merge, but this should not occur for 
                    ! a closed merge surface - throw error
                    call gdErrorHandler('MergeTMTubesTA: closed merge surface ' // & 
                        'detected, but all tubes are open. Not supported')

                end if 

            else
                ! Open-open merge
                call tmadaptor%MergeTMTubesOO(topomesh, hftubes, lftubes, mergefaces, &
                    wasmerged)

                ! Issue message in case of failure
                if (.not. wasmerged) then 
                    print *, 'MergeTMTubesTA: open-open merge not successful, continuing...'
                end if 

            end if 
        end if 

    end subroutine

    ! Closed-closed flux tube pair merging (non-separatrix)
    subroutine MergeTMTubesCCTA(tmadaptor, topomesh, tube1, tube2, mergefaces, &
        wasmerged)

        ! Description
        !============
        ! This routine performs the actual merge of two closed flux 
        ! tubes over their common faces given in 'mergefaces'. It is 
        ! assumed that all required checks on the input are done 
        ! beforehand in the calling function, in this case MergeTMTubes.
        ! See the algorithm section for more details. 

        ! Algorithm
        !==========
        ! 1) Determine all radial faces of both tubes
        ! 2) Check which radial faces are adjacent to one another. 
        ! 3) Radial faces that have an adjacent neighbour in the 
        ! other tube should be kept, others should be removed
        ! 4) Remove the merge faces and non-paired radial faces and
        ! rebuild the topomesh

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), intent(in)                 :: tube1, tube2
        integer(I8), dimension(:), intent(in)   :: mergefaces
        logical, intent(out)                    :: wasmerged

        ! Auxiliary
        logical, allocatable, dimension(:)      :: remf, remv
        integer(I8), allocatable, dimension(:)  :: tf1, tf2, tvf

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Initialize
        allocate(remf(topomesh%face%ntot), remv(topomesh%vert%ntot))
        remf = .false. 
        remv = .false.
        wasmerged = .false. 

        ! Get tube radial faces
        tf1 = topomesh%tube%GetFace(tube1)
        tf2 = topomesh%tube%GetFace(tube2)

        ! Mark faces for removal
        !=======================
        ! Mark merge faces for removal 
        remf(mergefaces) = .true. 
        
        ! Set removal to true for radial faces (will be set to false if pair found)
        remf(tf1) = .true.
        remf(tf2) = .true.

        ! Sanity check
        if (any(remf .and. topomesh%face%BF)) then
            ! Apparently some boundary faces were marked for deletion, 
            ! which should actually not happen. 
            call gdErrorHandler('MergeTMTubesCCTA: boundary faces were marked ' // &
                'for deletion, unexpected')
        end if 

        ! Check which radial faces to remove
        do i = 1, size(tf1)
            do j = 1, size(tf2)
                if (AreTMFacesAdjacent(topomesh%face, tf1(i), tf2(j))) then 
                    remf(tf1(i)) = .false. 
                    remf(tf2(j)) = .false.
                end if 
            end do 
        end do 

        ! Remove vertices that only have deleted faces
        do i = 1, topomesh%vert%ntot
            ! Get faces
            tvf = topomesh%vert%GetFace(i)

            if (all(remf(tvf))) then 
                remv(i) = .true.
            end if 
        end do

        ! Remove & rebuild
        !=================
        ! Remove faces
        call RemoveTopologicalMeshFaceLogical(topomesh, remf)
        call tmadaptor%RemoveFaceData(remf)

        ! Remove vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, remv)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Split faces if necessary
        call SplitTopologicalMeshFaces(topomesh, remf)  
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if  

        ! Recompute all interconnections, cells, etc
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Data (recompute)
        call AddTopologicalMeshData(topomesh)

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Remove faces from adaptor

        ! Set output
        wasmerged = .true.

    end subroutine

    ! Open-closed flux tube pair merging (non-separatrix)
    subroutine MergeTMTubesOCTA(tmadaptor, topomesh, tube1, tube2, mergefaces, &
        wasmerged)

        ! Description
        !============
        ! This routine performs the actual merge of an open and closed flux 
        ! tube over their common faces given in 'mergefaces'. It is 
        ! assumed that all required checks on the input are done 
        ! beforehand in the calling function, in this case MergeTMTubes.
        ! See the algorithm section for more details. 

        ! Algorithm
        !==========
        ! This is a very special case where the open tube either has
        ! two radial faces that connect in the same point with a radial
        ! face of the closed tube, or where these faces connect to a set 
        ! of aligned vessel parts. In either case, the radial faces of 
        ! the open tube are redefined to be aligned vessel parts and the
        ! common radial faces are removed. 

        ! 1) Determine which tube is open and take its radial faces
        ! 2) Redefine these faces as aligned vessel parts (and check if
        ! that's possible)
        ! 3) Remove all merge faces from the topomesh and any vertices
        ! that only connect to deleted faces
        ! 4) Rebuild the topomesh

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), intent(in)                 :: tube1, tube2
        integer(I8), dimension(:), intent(in)   :: mergefaces
        logical, intent(out)                    :: wasmerged

        ! Auxiliary
        logical, allocatable, dimension(:)      :: remf, remv
        integer(I8), allocatable, dimension(:)  :: tf, tvf

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Initialize
        wasmerged = .false. 
        allocate(remf(topomesh%face%ntot), remv(topomesh%vert%ntot))
        remf = .false. 
        remv = .false.

        ! Get open tube radial faces
        if (topomesh%tube%isclosed(tube1)) then 
            tf = topomesh%tube%GetFace(tube2)
        elseif (topomesh%tube%isclosed(tube2)) then 
            tf = topomesh%tube%GetFace(tube1)
        else
            ! Shouldn't happen
            call gdErrorHandler('MergeTMTubesOC: both tubes are open, ' // & 
                'unexpected')
        end if

        ! Mark faces for removal
        !=======================
        ! Retype open tube radial faces
        if (.not. all(topomesh%face%type(tf([1, size(tf)])) == TMfacebndID)) then 
            ! We can only retype to aligned vessel parts if they're actually
            ! vessel parts. If not, throw error
            call gdErrorHandler('MergeTMTubesOC: not all radial faces of' // & 
                'open tube are vessel boundary faces, unexpected')
        end if 
        topomesh%face%type(tf([1, size(tf)])) = TMfacealbndID 

        ! Assign flux surface IDs (just take one of the vertices...)
        topomesh%face%fsID(tf([1, size(tf)])) = topomesh%vert%fsID(topomesh%face%vert(tf(1), 1))

        ! Mark merge faces for removal
        remf(mergefaces) = .true. 
        where (topomesh%face%BF) remf = .false. 
        
        ! Remove vertices that only have deleted faces
        do i = 1, topomesh%vert%ntot
            ! Get faces
            tvf = topomesh%vert%GetFace(i)

            if (all(remf(tvf))) then 
                remv(i) = .true.
            end if 
        end do

        ! Remove & rebuild
        !=================
        ! Remove faces
        call RemoveTopologicalMeshFaceLogical(topomesh, remf)
        call tmadaptor%RemoveFaceData(remf)

        ! Remove vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, remv)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Split faces if necessary
        call SplitTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if    

        ! Recompute all interconnections, cells, etc
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Data (recompute)
        call AddTopologicalMeshData(topomesh)

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Remove faces from adaptor

        ! Set output
        wasmerged = .true.

    end subroutine

    ! Open-open flux tube pair merging (non-separatrix)
    subroutine MergeTMTubesOOTA(tmadaptor, topomesh, tube1, tube2, mergefaces, &
        wasmerged)

        ! Description
        !============
        ! This routine performs the actual merge of multiple open flux 
        ! tubes over their common faces given in 'mergefaces'. It is 
        ! assumed that all required checks on the input are done 
        ! beforehand in the calling function, in this case MergeTMTubes.
        ! See the algorithm section for more details.
        
        ! Algorithm
        !==========
        ! Here, the merge faces should form a simple, open polygon. 
        ! Furthermore, it is assumed that these faces are sorted from 
        ! start to end, such that we can determine the starting and 
        ! ending vertices, as well as the adjacent radial faces. Basically, 
        ! we try to keep as much of the vessel boundary structure to be 
        ! defined as vessel and only redefine vessel parts to be aligned
        ! if it is not otherwise possible. In practice, this means that 
        ! any vessel boundary that does not connect only through vessel boundaries
        ! to the outer vertices of the merging faces has to be redefined
        ! to an aligned vessel type. Otherwise, the merge will not result 
        ! in a conforming flux tube. 

        ! 1) Determine all faces that will be part of the final merged 
        ! tube (these are basically all the non-merging faces and any 
        ! radial faces at outer sides of each tube. 'Inner' radial faces
        ! will be removed, if they are even present). These faces should
        ! form a closed polygon.
        ! 1) Determine the outer vertices of the merge faces
        ! 2) Start at one outer vertex and walk along the faces of the 
        ! closed tube until we either reach an aligned part or the next 
        ! vertex. In the former case, all intermediate vertices should be 
        ! marked for removal. In the latter case, at least one tangency
        ! point should be retained (keep the one that gives the largest
        ! difference in psi values between the merging face psi values).
        ! Repeat this procedure for the second outer vertex. Mark which 
        ! vessel faces should be kept as vessel faces.
        ! 3) Retype non-marked vessel faces as aligned vessel parts
        ! 4) Remove merge faces and vertices with only removed faces, 
        ! but don't remove any boundary merge faces
        ! 5) Rebuild the topomesh

        ! Notes
        !======
        ! Note 1: in this type of merge, it is actually possible to 
        ! form tubes that have overlapping psi values between high and
        ! low flux side. This is not desireable, since these tubes 
        ! have maximally zero radial length/delta psi and will therefore
        ! always pop up in a merge, leading to excessive tube merging.
        ! Hence, it also leads to excessively coarse grids. Therefore, 
        ! we will not perform the merge if this case is detected 
        ! (a message will be shown). Flux surface IDs of the merge 
        ! faces will be added to the illegal flux surfaces to merge 
        ! over (these should not change during topomesh adaptations)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), dimension(:), intent(in)   :: mergefaces, tube1, &
            tube2
        logical, intent(out)                    :: wasmerged

        ! Auxiliary
        logical, allocatable, dimension(:)      :: remf, remv
        integer(I8)                             :: ne, ind1, ind2, &
            thisf, v1, v2, si, ei
        integer(I8), allocatable, dimension(:)  :: tf, tfb, &
            tempi, tfbv, mergevert, sortind, &
            polygonID, bndvert1, bndvert2, bndvert3, bndvert4, &
            tvf, tempf, avpfsID, tv, tvfsID, commonmergefaces
        integer(I8), allocatable, dimension(:, :)   :: tfv
        logical                                 :: noalignedfacev1tov2, &
            noalignedfacev2tov1
        logical, allocatable, dimension(:)      :: ispolygonstart, &
            isbranchingpolygon, isavp, retypevert, isbndf1, isbndf2
        real(R8), allocatable, dimension(:)     :: mergepsival, tpsi, &
            dpsi, tvfval
        type(IntegerDynamicArrayUDT)            :: tfbida, tfida1, tfida2

        ! Loop
        integer(I8)                             :: i, k

        ! Initialize
        !===========
        ! Initialize
        allocate(remf(topomesh%face%ntot), remv(topomesh%vert%ntot))
        remf = .false. 
        remv = .false.
        mergepsival = topomesh%fsfval%Get(topomesh%face%fsID(mergefaces))
        wasmerged = .false. 

        ! Construct closed polygon
        !=========================
        ! For each tube, add bounding faces 
        tfbida = ConstructIntegerDynamicArray()
        tfida1 = ConstructIntegerDynamicArray()
        tfida2 = ConstructIntegerDynamicArray()
        do i = 1, size(tube1)
            ! Radial faces: only first and last
            tf = topomesh%tube%GetFace(tube1(i))
            call tfbida%Append(tf([1, size(tf)]))

            ! Boundary faces: currently all (non-boundary merge faces will be deleted afterwards)
            tf = topomesh%tube%GetBndFace(tube1(i), 1)
            call tfbida%Append(tf)
            call tfida1%Append(tf)
            tf = topomesh%tube%GetBndFace(tube1(i), 2)
            call tfbida%Append(tf)
            call tfida1%Append(tf)
        end do 
        do i = 1, size(tube2)
            ! Radial faces: only first and last
            tf = topomesh%tube%GetFace(tube2(i))
            call tfbida%Append(tf([1, size(tf)]))

            ! Boundary faces: currently all (non-boundary merge faces will be deleted afterwards)
            tf = topomesh%tube%GetBndFace(tube2(i), 1)
            call tfbida%Append(tf)
            call tfida2%Append(tf)
            tf = topomesh%tube%GetBndFace(tube2(i), 2)
            call tfbida%Append(tf)
            call tfida2%Append(tf)
        end do 

        ! Remove non-common merge faces
        commonmergefaces = GetCommonElements(&
            GetCommonElements(tfida1%Get(), mergefaces), &
            GetCommonElements(tfida2%Get(), mergefaces))

        tempi = tfbida%Get()
        !allocate(nonbndmergefaces(count(.not. topomesh%face%BF(mergefaces))))
        !nonbndmergefaces = pack(mergefaces, .not. topomesh%face%BF(mergefaces))
        !call SetDiff(tempi, nonbndmergefaces, tfb)
        call SetDiff(tempi, commonmergefaces, tfb)

        ! Sort faces
        tfv = topomesh%face%vert(tfb, :)
        ne = size(tfv, 1)
        allocate(sortind(ne), ispolygonstart(ne), isbranchingpolygon(ne))
        call SortPolygonEdges(tfv, ne, sortind, &
            ispolygonstart, isbranchingpolygon, polygonID)
        tfv(:, 1) = tfv(sortind, 1)
        tfv(:, 2) = tfv(sortind, 2)
        tfb = tfb(sortind)
        deallocate(sortind)

        ! Check if single simple polygon
        if ((maxval(polygonID)) > 1 ) then 
            call gdErrorHandler('MergeTMTubesOOTA: tube polygon consists of ' // &
                'multiple polygons. Unexpected and unsupported')
        elseif (any(isbranchingpolygon)) then 
            ! Try with dedicated method
            tempf = GetClosedPolygonFromTopomeshFaces(topomesh, tfb)
            if (size(tempf) == 0) then 
                call WriteTopologicalMesh(topomesh, 'topomesh_error.dat')
                call gdErrorHandler('MergeTMTubesOOTA: tube polygon likely ' // & 
                    'contains separatrix part, but could not form closed ' // & 
                    'polygon, even with dedicated routine. Unexpected, check input')
            end if 
            tfb = tempf
            tfv = topomesh%face%vert(tfb, :)
        end if 

        ! Extract vertices
        call ExtractPolygonVertices(topomesh%face%vert(tfb, :), &
            size(tfb), tfbv)

        ! Check if closed polygon
        if (tfbv(1) /= tfbv(size(tfbv))) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
            call gdErrorHandler('MergeTMTubesOOTA: tube polygon is not ' // & 
                'closed, unexpected')
        end if 

        ! Determine outer vertices
        !=========================
        ! Extract vertices (faces assumed sorted before)
        call ExtractPolygonVertices(topomesh%face%vert(mergefaces, :), &
            size(mergefaces), mergevert)

        ! Sanity check
        if (mergevert(1) == mergevert(size(mergevert))) then 
            call gdErrorHandler('MergeTMTubesOOTA: merge faces seem to ' // & 
                'form closed polygon, unexpected here')
        end if 

        ! Determine initial outer vertices
        v1 = mergevert(1)
        v2 = mergevert(size(mergevert))

        ! Check if we need to move the vertices further along the merge
        ! faces

        ! Determine aligned parts
        !========================
        ! Initialize
        noalignedfacev1tov2 = .false.
        noalignedfacev2tov1 = .false.
        allocate(bndvert1(0), bndvert2(0), bndvert3(0), bndvert4(0))
        allocate(isavp(size(tfb)), avpfsID(size(tfb)), isbndf1(size(tfb)), &
            isbndf2(size(tfb)))
        isavp = .true. ! will be set to false where applicable later
        isbndf1 = .true.
        isbndf2 = .true. 
        avpfsID = 0 ! flux surface ID for aligned vessel parts
        where (.not. (topomesh%face%type(tfb) == TMfacebndID)) isavp = .false. ! also ignore non-vessel boundaries

        ! Find v1, v2 in closed polygon 
        ind1 = findloc(tfbv, v1, 1, back=.false.)
        ind2 = findloc(tfbv, v2, 1, back=.false.)
        if (ind1 == 0) then 
            call gdErrorHandler('MergeTMTubesOO: first outer vertex does ' // & 
                'not appear in tube polygon, unexpected')
        end if
        if (ind2 == 0) then 
            call gdErrorHandler('MergeTMTubesOO: second outer vertex does ' // & 
                'not appear in tube polygon, unexpected')
        end if

        ! Set logicals determining boundary faces to false where necessary
        if (ind1 < ind2) then 
            isbndf1(1:max(ind1-1, 1)) = .false.
            isbndf1(min(ind2, size(tfb)):) = .false.
            isbndf2 = .not. isbndf1
        else
            isbndf2(1:max(ind2-1, 1)) = .false.
            isbndf2(min(ind1, size(tfb)):) = .false.
            isbndf1= .not. isbndf2
        end if 

        ! Step 'forward' for v1 
        k = ind1
        bndvert1 = [bndvert1, tfbv(k)]
        do while (.true.)
            ! Take the current face
            if (k == size(tfbv)) then 
                k = 1
            end if 
            thisf = tfb(k)

            ! Check if it is a vessel boundary
            if (topomesh%face%type(thisf) == TMfacebndID) then 
                isavp(k) = .false.
                isbndf1(k) = .false.  
            elseif (any(topomesh%face%type(thisf) == TMfacealignedID)) then 
                ! Aligned boundary found, exit
                exit
            end if 

            ! Update k 
            k = k + 1

            ! Add vertex
            bndvert1 = [bndvert1, tfbv(k)]

            ! Check if we encountered v2
            if (tfbv(k) == v2) then 
                ! No aligned face found between v1 and v2 
                noalignedfacev1tov2 = .true.
                exit
            end if 
        end do 

        ! Step 'backward' for v1 
        k = ind1
        bndvert2 = [bndvert2, tfbv(k)]
        do while (.true.)
            ! Update k 
            k = k - 1

            ! Take the current face
            if (k == 0) then 
                k = size(tfbv) - 1
            end if 
            thisf = tfb(k)

            ! Check if it is a vessel boundary
            if (topomesh%face%type(thisf) == TMfacebndID) then 
                isavp(k) = .false. 
                isbndf2(k) = .false. ! stepping along other side now!
            elseif (any(topomesh%face%type(thisf) == TMfacealignedID)) then 
                ! Aligned boundary found, exit
                exit
            end if 

            ! Add vertex
            bndvert2 = [bndvert2, tfbv(k)]

            ! Check if we encountered v2
            if (tfbv(k) == v2) then 
                ! No aligned face found between v1 and v2 
                noalignedfacev2tov1 = .true.
                exit
            end if 
        end do 

        ! Step 'forward' for v2
        k = ind2
        bndvert3 = [bndvert3, tfbv(k)]
        do while (.true. .and. .not. (noalignedfacev2tov1))
            ! Take the current face
            if (k == size(tfbv)) then 
                k = 1
            end if 
            thisf = tfb(k)

            ! Check if it is a vessel boundary
            if (topomesh%face%type(thisf) == TMfacebndID) then 
                isavp(k) = .false. 
                isbndf2(k) = .false. 
            elseif (any(topomesh%face%type(thisf) == TMfacealignedID)) then 
                ! Aligned boundary found, exit
                exit
            end if 

            ! Update k 
            k = k + 1

            ! Add vertex
            bndvert3 = [bndvert3, tfbv(k)]

            ! Check if we encountered v1
            if (tfbv(k) == v1) then 
                ! No aligned face found between v1 and v2 
                noalignedfacev2tov1 = .true.
                exit
            end if 
        end do 

        ! Step 'backward' for v2 
        k = ind2
        bndvert4 = [bndvert4, tfbv(k)]
        do while (.true. .and. .not. (noalignedfacev1tov2))
            ! Update k 
            k = k - 1

            ! Take the current face
            if (k == 0) then 
                k = size(tfbv) - 1
            end if 
            thisf = tfb(k)

            ! Check if it is a vessel boundary
            if (topomesh%face%type(thisf) == TMfacebndID) then 
                isavp(k) = .false. 
                isbndf1(k) = .false. ! Stepping along other side now!
            elseif (any(topomesh%face%type(thisf) == TMfacealignedID)) then 
                ! Aligned boundary found, exit
                exit
            end if 

            ! Add vertex
            bndvert4 = [bndvert4, tfbv(k)]

            ! Check if we encountered v1
            if (tfbv(k) == v1) then 
                ! No aligned face found between v1 and v2 
                noalignedfacev1tov2 = .true.
                exit
            end if 
        end do 

        ! Sanity check
        if (noalignedfacev1tov2 .and. noalignedfacev2tov1) then 
            ! This indicates a flux tube with only vessel faces, which is
            ! currently not supported (not really an issue for this 
            ! routine but it is for grid generation later on)
            call gdErrorHandler('MergeTMTubesOOTA: tube merging would result ' // & 
                'in tube with no aligned faces, not yet supported')
        end if  

        ! Retype faces
        !=============
        ! Check which vessel parts to align and determine flux surface IDs
        if (any(isavp)) then 
            si = 0
            ei = 0
            do while (.true.)
                ! Find avp segment
                si = findloc(isavp(ei+1:), .true., 1, back=.false.) + ei
                if (si == ei) then 
                    ! No more parts found, exit
                    exit 
                end if 
                ei = findloc(isavp(si:), .false., 1, back=.false.) + si - 2
                if (ei == si-2) then 
                    ei = size(isavp)
                end if 

                ! Set type
                topomesh%face%type(tfb(si:ei)) = TMfacealbndID

                ! Determine flux surface ID based on flux values of vertices
                call SetDiff([topomesh%face%vert(tfb(si:ei), 1), topomesh%face%vert(tfb(si:ei), 2)], &
                    [topomesh%face%vert(commonmergefaces, 1), topomesh%face%vert(commonmergefaces, 2)], tv)
                tvfsID = topomesh%vert%fsID(tv)

                ! Keep only vertices with flux surface ID
                tv = pack(tv, tvfsID /= 0)
                tvfsID = pack(tvfsID, tvfsID /= 0)
                tvfval = topomesh%fsfval%Get(tvfsID)

                if (size(tvfval) == 0) then 
                    call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                    print *, 'vertices: ', tv
                end if 

                
                ! Compute minimal distance in terms of psi value w.r.t.
                ! the merging faces & set flux surface ID
                if (abs(maxval(tvfval) - minval(mergepsival)) <  abs(minval(tvfval) - maxval(mergepsival))) then 
                    topomesh%face%fsID(tfb(si:ei)) = tvfsID(maxloc(tvfval, 1))
                else
                    topomesh%face%fsID(tfb(si:ei)) = tvfsID(minloc(tvfval, 1))
                end if
                
            end do 
        end if 
        
        ! Retype vertices 
        !================
        ! First part
        if (noalignedfacev1tov2) then 
            ! Need to keep one tangency point of type 1 at this side, 
            ! take one with maximal psi difference

            ! Sanity check
            if (all(topomesh%vert%type(bndvert1) /= TMvertextp1ID)) then 
                call WriteTopologicalMesh(topomesh, 'topomesh_error')
                print *, 'vertices: ', bndvert1
                call gdErrorHandler('MergeTMTubesOO: expected to have ' // & 
                    'at least one tangency point type 1 in boundary but ' // & 
                    'found none')
            end if 

            ! Determine tangency point type 1 vertex with maximal distance
            ! to merge boundary
            tpsi = topomesh%vert%fval(bndvert1)
            dpsi = abs(maxval(mergepsival) - tpsi)
            where (topomesh%vert%type(bndvert1) /= TMvertextp1ID) dpsi = -posinfval_R8()
            
            ! Retype
            allocate(retypevert(size(bndvert1)))
            retypevert = .true. !topomesh%vert%type(bndvert1) == TMvertextp1ID ! retype all for now
            retypevert(maxloc(dpsi)) = .false. ! exclude tangency point
            where (retypevert) topomesh%vert%type(bndvert1) = TMvertexbndID
            deallocate(retypevert) 

        else
            ! Need to retype all 'inner' boundary vertices to regular 
            ! boundary vertices for removal later on

            ! Bndvert 1, 4 (does not go to v2, so we can retype here)
            topomesh%vert%type(bndvert1(2:size(bndvert1)-1)) = TMvertexbndID
            topomesh%vert%type(bndvert4(2:size(bndvert4)-1)) = TMvertexbndID
        end if 

        ! Second part
        if (noalignedfacev2tov1) then 
            ! Need to keep one tangency point of type 1 at this side, 
            ! take one with maximal psi difference

            ! Sanity check
            if (all(topomesh%vert%type(bndvert2) /= TMvertextp1ID)) then 
                call WriteTopologicalMesh(topomesh, 'topomesh_error')
                print *, 'vertices: ', bndvert2
                call gdErrorHandler('MergeTMTubesOO: expected to have ' // & 
                    'at least one tangency point type 1 in boundary but ' // & 
                    'found none')
            end if 

            ! Determine tangency point type 1 vertex with maximal distance
            ! to merge boundary
            tpsi = topomesh%vert%fval(bndvert2)
            dpsi = abs(maxval(mergepsival) - tpsi)
            where (topomesh%vert%type(bndvert2) /= TMvertextp1ID) dpsi = -posinfval_R8()
            
            ! Retype
            allocate(retypevert(size(bndvert2)))
            retypevert = .true. !topomesh%vert%type(bndvert2) == TMvertextp1ID ! retype all for now
            retypevert(maxloc(dpsi)) = .false. ! Exclude tangency point
            where (retypevert) topomesh%vert%type(bndvert2) = TMvertexbndID
            deallocate(retypevert) 


        else
            ! Need to retype all 'inner' boundary vertices to regular 
            ! boundary vertices for removal later on 

            ! Bndvert 2, 3 (does not go to v2, so we can retype here)
            topomesh%vert%type(bndvert2(2:size(bndvert2)-1)) = TMvertexbndID
            topomesh%vert%type(bndvert3(2:size(bndvert3)-1)) = TMvertexbndID
        end if 

        ! Mark faces for removal
        !=======================
        ! Mark merge faces for removal that are not boundary faces and 
        ! that are common between both tubes
        where (.not. topomesh%face%BF(commonmergefaces)) remf(commonmergefaces) = .true.
        
        ! Mark any non-boundary radial faces of tubes to be removed
        do i = 1, size(tube1)
            tf = topomesh%tube%GetFace(tube1(i))
            remf(tf(2:size(tf)-1)) = .true.
        end do 
        do i = 1, size(tube2)
            tf = topomesh%tube%GetFace(tube2(i))
            remf(tf(2:size(tf)-1)) = .true.
        end do 

        ! Remove vertices that only have deleted faces
        do i = 1, topomesh%vert%ntot
            ! Get faces
            tvf = topomesh%vert%GetFace(i)

            if (all(remf(tvf))) then 
                remv(i) = .true.
            end if 
        end do

        ! Remove & rebuild
        !=================
        ! Remove faces
        call RemoveTopologicalMeshFaceLogical(topomesh, remf)
        call tmadaptor%RemoveFaceDataLogical(remf)

        ! Remove vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, remv)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Split faces if necessary
        call SplitTopologicalMeshFaces(topomesh, remf)   
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 
        ! Recompute all interconnections, cells, etc
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Data (recompute)
        call AddTopologicalMeshData(topomesh)

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Adjust adaptor
        

        ! Set output
        wasmerged = .true.

    end subroutine

    ! Separatrix flux tube pair merging 
    subroutine MergeTMTubesSTA(tmadaptor, topomesh, tube1, tube2, &
        mergefaces, wasmerged)

        ! Description
        !============
        ! This routine performs the actual merge over a separatrix boundary
        ! of which the faces are given in 'mergefaces'. It is 
        ! assumed that all required checks on the input are done 
        ! beforehand in the calling function, in this case MergeTMTubes.
        ! The merge algorithm(s) presented here are rather involved since
        ! merging away separatrices is, at best, difficult. Therefore, 
        ! there is a non-negligible possibility that this routine will 
        ! fail and not merge away the tubes. Currently, not all 
        ! separatrix merging cases are supported. The ones that are, 
        ! are the following:
        ! - merging where at least on one side all tubes have separatrices.
        !   this is called a 'collapse' since we're basically collapsing
        !   that tube. This is a very useful merge operation to simplify
        !   topologies with x-points that lie very close together. In 
        !   principle, this can only happen for a single tube. 

        ! Algorithm
        !==========
        ! Collapse merge
        !---------------
        ! This merge can only proceed if one tube set is fully contained
        ! by separatrix segments, and if that tube set has at least one
        ! x-point in the low and high flux boundary. Essentially, both 
        ! separatrices are merged together by taking pieces of both. All
        ! x-points are retained. See the dedicated subroutine for 
        ! more information

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), dimension(:), intent(in)   :: mergefaces, tube1, tube2
        logical                                 :: wasmerged

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: hfbndf, lfbndf, &
            hfbndv, lfbndv, tubemergefaces
        logical                                 :: iscollapsemerge, &
            ishfxpfound, islfxpfound, isallsep
        logical, allocatable, dimension(:)      :: collapsetube1, &
            collapsetube2, ismergeface

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Set output
        wasmerged = .false. 

        ! Associate for ease
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            tube        => topomesh%tube    &
            )

        ! Determine possible merge cases
        !===============================
        ! Collapse merge
        !---------------
        ! Initialize switch
        iscollapsemerge = .false.

        ! Initialize additional data
        collapsetube1 = spread(.false., 1, size(tube1)) 
        collapsetube2 = spread(.false., 1, size(tube2))

        ! Check tube1
        do i = 1, size(tube1)
            ! Initialize
            ishfxpfound = .false.
            islfxpfound = .false. 
            isallsep = .true.  

            ! Get boundary faces & vertices
            hfbndf = tube%GetHighFluxBndFace(tube1(i))
            lfbndf = tube%GetLowFluxBndFace(tube1(i))
            hfbndv = [face%vert(hfbndf, 1), face%vert(hfbndf, 2)]
            lfbndv = [face%vert(lfbndf, 1), face%vert(lfbndf, 2)]

            ! Check if all boundaries are separatrices
            if (.not. all(face%type([hfbndf, lfbndf]) == TMfacesepID)) then 
                isallsep = .false.
            end if 

            ! Check if an X-point is found at high and low field side
            if (any(vert%type(hfbndv) == TMvertexsaddleID)) then 
                ishfxpfound = .true.
            end if 
            if (any(vert%type(lfbndv) == TMvertexsaddleID)) then 
                islfxpfound = .true.
            end if
            
            ! Check
            if (ishfxpfound .and. islfxpfound .and. isallsep) then 
                collapsetube1(i) = .true.
            end if
        end do 
        

        ! Check tube2
        do i = 1, size(tube2)
            ! Initialize
            ishfxpfound = .false.
            islfxpfound = .false. 
            isallsep = .true. 

            ! Get boundary faces & vertices
            hfbndf = tube%GetHighFluxBndFace(tube2(i))
            lfbndf = tube%GetLowFluxBndFace(tube2(i))
            hfbndv = [face%vert(hfbndf, 1), face%vert(hfbndf, 2)]
            lfbndv = [face%vert(lfbndf, 1), face%vert(lfbndf, 2)]

            ! Check if all boundaries are separatrices
            if (.not. all(face%type([hfbndf, lfbndf]) == TMfacesepID)) then 
                isallsep = .false.
                exit 
            end if 

            ! Check if an X-point is found at high and low field side
            if (any(vert%type(hfbndv) == TMvertexsaddleID)) then 
                ishfxpfound = .true.
            end if 
            if (any(vert%type(lfbndv) == TMvertexsaddleID)) then 
                islfxpfound = .true.
            end if 

            ! Check
            if (ishfxpfound .and. islfxpfound .and. isallsep) then 
                collapsetube2(i) = .true.
            end if
        end do 
        

        ! Check final conditions
        if (any(collapsetube1) .or. any(collapsetube2)) then 
            iscollapsemerge = .true.
        end if 

        ! Housekeeping
        !=============
        end associate

        ! Merge
        !======
        if (iscollapsemerge) then 
            ! Loop until either a tube was merged or until all tubes
            ! were considered
            allocate(ismergeface(topomesh%face%ntot))
            ismergeface = .false.
            ismergeface(mergefaces) = .true. 
            do i = 1, size(tube1)
                if (collapsetube1(i)) then 
                    ! Determine the current merge faces
                    hfbndf = topomesh%tube%GetHighFluxBndFace(tube1(i))
                    lfbndf = topomesh%tube%GetLowFluxBndFace(tube1(i))
                    if (all(ismergeface(hfbndf)) .and. .not. any(ismergeface(lfbndf))) then 
                        tubemergefaces = hfbndf 
                    elseif (all(ismergeface(lfbndf)) .and. .not. any(ismergeface(hfbndf))) then
                        tubemergefaces = lfbndf
                    else
                        ! Seems like neither side is fully marked as merge faces, unexpected. Throw error
                        call gdErrorHandler('MergeTMTubesSTA: neither ' // & 
                            'high nor low flux side tube faces are all marked ' // & 
                            'as merging faces, unexpected')
                    end if  

                    ! Merge tube
                    call tmadaptor%CollapseTMSeparatrixTube(topomesh, &
                        tube1(i), tubemergefaces, wasmerged)

                    ! Check if we merged, if so, exit
                    if (wasmerged) then 
                        exit 
                    end if 
                end if 
            end do
            if (.not. wasmerged) then 
                do i = 1, size(tube2)
                    if (collapsetube2(i)) then 
                        ! Determine the current merge faces
                        hfbndf = topomesh%tube%GetHighFluxBndFace(tube2(i))
                        lfbndf = topomesh%tube%GetLowFluxBndFace(tube2(i))
                        if (all(ismergeface(hfbndf)) .and. .not. any(ismergeface(lfbndf))) then 
                            tubemergefaces = hfbndf 
                        elseif (all(ismergeface(lfbndf)) .and. .not. any(ismergeface(hfbndf))) then
                            tubemergefaces = lfbndf
                        else
                            ! Seems like neither side is fully marked as merge faces, unexpected. Throw error
                            call gdErrorHandler('MergeTMTubesSTA: neither ' // & 
                                'high nor low flux side tube faces are all marked ' // & 
                                'as merging faces, unexpected')
                        end if  

                        ! Merge tube
                        call tmadaptor%CollapseTMSeparatrixTube(topomesh, &
                            tube2(i), tubemergefaces, wasmerged)

                        ! Check if we merged, if so, exit
                        if (wasmerged) then 
                            exit 
                        end if 
                    end if 
                end do 
            end if 
        end if

    end subroutine

    ! Separatrix tube collapsing
    subroutine CollapseTMSeparatrixTubeTA(tmadaptor, topomesh, mergetube, &
        mergefaces, wasmerged)

        ! Description
        !============
        ! This routine 'collapses' a tube that has only separatrix
        ! boundaries and at least one x-point on both hf and lf side. It
        ! is assumed that these conditions have been checked beforehand 
        ! (this routine should only be called from MergeTMTubesSTA). If
        ! this is not the case, the merge is likely not going to be 
        ! successful. The tube is collapsed by removing all faces in 
        ! mergefaces eventually while retaining the x-point. The latter
        ! has to be inserted in some way in the other separatrix segment. 
        ! This may cause quite some local distortion depending on how 
        ! far away this x-point is from the other separatrix. In some
        ! (exceptional) cases this insertion may fail, which will abort
        ! the merge. 

        ! Algorithm
        !==========
        ! 1)    Extract the non-merged separatrix from the tube
        ! 2)    Find the x-point(s) to be inserted in the non-merged 
        !       separatrix segment
        ! 3)    Insert the x-point(s) by finding the intersection between
        !       its streamline(s) and the non-merged separatrix. If the 
        !       same x-point appears multiple times in the separatrix, 
        !       then it has to be inserted as many times as it appears
        !       (e.g. when the x-point is close to a core region). If
        !       the streamline intersects with a vessel boundary of the
        !       tube, abort the merge.
        ! 4)    If desired/possible/necessary, smooth/simplify the 
        !       separatrix near the inserted x-point(s)
        ! 5)    Remove the merge faces and rebuild the topomesh

        ! Note: in step 3, there may be alternatives in order to let the
        ! merge proceed. For example, one could insert an aligned vessel
        ! part between the intersection and the next separatrix and
        ! consider that to be a separatrix segment. This is not yet 
        ! implemented as it is quite unlikely for current configurations
        ! to appear (although it was also quite unlikely a few weeks ago
        ! that this routine would be necessary, so, yeah.)

        ! Note: currently, no actual streamlines are traced, but a 
        ! line segment between the x-point and the closes point of the 
        ! separatrix that will be retained is determined. We could
        ! also trace streamlines, but typically these are quite badly
        ! behaved and hard to handle near X-points. 

        ! Note: we explicitly do *not* update the separatrix face IDs to
        ! be consistent after merging, as this may lead to issues 
        ! when determining the high and low field side of the tube. This
        ! is also the most consistent approach, since the flux surface 
        ! IDs are, in principle, not the same and also the values 
        ! are different. 


        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor 
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), dimension(:), intent(in)   :: mergefaces
        integer(I8), intent(in)                 :: mergetube
        logical, intent(out)                    :: wasmerged

        ! Auxiliary
        integer(I8)                             :: tracedir, tf1, tf2, &
            nint, nsl, ntf, nxp, newfsID, ind1, ind2
        integer(I8), allocatable, dimension(:)  :: allhf, alllf, &
            sepfaces, thf, tlf, mergevert, mergesepfaces, tubefaces, &
            radialbndf, s1, s2, sepfaceID, xpointID, news2, usepfaceID, &
            txp, tnews2, sortindex, v1, v2, tvf, tvfna, tuberadface
        logical                                 :: ishfsep, islfsep, &
            isintersectionfound
        logical, allocatable, dimension(:)      :: ismergeface, &
            isseparatrixface, ishfface, islfface, remf, remv, &
            istuberadface
        real(R8)                                :: x1, y1, x2, y2,&
            xb(1:2), yb(1:2), fsfval, xxp, yxp, mindist, tdist, tx1, ty1, &
            tx2, ty2
        real(R8), allocatable, dimension(:)     :: xt, yt, xs, ys, xf, &
            yf, xint, yint, s1r, s2r, newx, newy, news2r, tnews2r, ds, &
            txv, tyv
        type(StreamlineUDT), allocatable        :: streamlines(:)
        type(RealDynamicArrayUDT), allocatable      :: xintrda(:, :), yintrda(:, :), &
            s1rrda(:, :), s2rrda(:, :), xfnew(:), yfnew(:)
        type(IntegerDynamicArrayUDT), allocatable   :: s1ida(:, :), &
            s2ida(:, :)

        ! Loop
        integer(I8)                             :: i, j, cc

        ! Initialize
        !===========
        ! Set output
        wasmerged = .false. 

        ! Sanity checks
        if (mergetube == 0) then 
            print *, 'CollapseTMSeparatrixTubeTA: no tubes given as input, returning...'
            wasmerged = .true.
            return 
        end if 
        if (size(mergefaces) == 0) then 
            call gdErrorHandler('CollapseTMSeparatrixTubeTA: no merge faces ' // & 
                'given, cannot proceed')
        end if

        ! Unpack for ease
        associate(&
            vert    => topomesh%vert,   &
            face    => topomesh%face,   &
            tube    => topomesh%tube    &
            )

        ! Set some useful logicals
        allocate(ismergeface(face%ntot), isseparatrixface(face%ntot), &
            istuberadface(face%ntot))
        ismergeface = .false.
        isseparatrixface = .false. 
        istuberadface = .false. 
        ismergeface(mergefaces) = .true.
        tuberadface = tube%GetFace(mergetube)
        istuberadface(tuberadface) = .true. 

        ! Extract separatrix 
        !===================
        ! Check high flux side
        allocate(allhf(0))
        ishfsep = .true.
        thf = tube%GetHighFluxBndFace(mergetube)
        allhf = [allhf, thf]

        ! Check
        if (any(ismergeface(thf))) then 
            ishfsep = .false. 
        end if 

        ! Check low flux side
        allocate(alllf(0))
        islfsep = .true.
        tlf = tube%GetLowFluxBndFace(mergetube)
        alllf = [alllf, tlf]

        ! Check
        if (any(ismergeface(tlf))) then 
            islfsep = .false. 
        end if 

        ! Sanity checks
        if (ishfsep .and. islfsep) then 
            ! This indicates something wrong in the input
            call gdErrorHandler('CollapseTMSeparatrixTubeTA: both high and ' // & 
                'low flux faces have been identified as final separatrix, unexpected')
        elseif ((.not. ishfsep) .and. (.not. islfsep)) then 
            ! This indicates something wrong in the input
            call gdErrorHandler('CollapseTMSeparatrixTubeTA: neither high and ' // & 
                'low flux faces have been identified as final separatrix, unexpected')
        end if 
        if (ishfsep) then 
            if (.not. all(ismergeface(alllf))) then 
                ! Normally all low flux surfaces should be the merge 
                ! faces now
                call gdErrorHandler('CollapseTMSeparatrixTubeTA: not all ' // & 
                    'low flux faces are merge faces though high flux faces ' // & 
                    'are non-merging separatrix faces, unexpected')
            end if
        end if
        if (islfsep) then 
            if (.not. all(ismergeface(allhf))) then 
                ! Normally all high flux surfaces should be the merge 
                ! faces now
                call gdErrorHandler('CollapseTMSeparatrixTubeTA: not all ' // & 
                    'high flux faces are merge faces though low flux faces ' // & 
                    'are non-merging separatrix faces, unexpected')
            end if
        end if

        ! Extract vessel faces (has to be two or zero)
        radialbndf = pack(tuberadface, face%type(tuberadface) == TMfacebndID)

        ! Sanity checks
        if (size(radialbndf) == 0) then 
            ! The tubes form a closed region  - take any radial face
            radialbndf = [radialbndf, tuberadface(1)]
        elseif (size(radialbndf) == 2) then 
            ! all good
        else
            ! This indicates some complicated topology which we don't yet
            ! support
            print *, 'CollapseTMSeparatrixTubeTA: could not exactly find ' // & 
                'zero or two vessel boundary faces, unexpected. Aborting merge'
            return 
        end if 

        ! Orden the faces of the tube
        tubefaces = GetClosedPolygonFromTopomeshFaces(topomesh, [radialbndf, allhf, alllf])

        ! Extract (sorted) separatrix faces - need to flip one of the two
        ! to have them oriented in the same direction
        allocate(islfface(face%ntot), ishfface(face%ntot))
        islfface = .false.
        ishfface = .false.
        islfface(alllf) = .true.
        ishfface(allhf) = .true.
        
        ! Sanity checks
        if (all(ishfface(tubefaces)) .or. all(islfface(tubefaces))) then 
            call gdErrorHandler('CollapseSeparatrixTubeTA: all faces are ' // & 
                'either high or low flux faces, unexpected.')
        end if 
        if (ishfsep) then 
            ! Hedge for faces not being stored continuously in mergefaces
            if (islfface(tubefaces(1)) .and. islfface(tubefaces(size(tubefaces)))) then 
                ind1 = findloc(.not. islfface(tubefaces), .true., 1, back=.false.)-1
                ind2 = findloc(.not. islfface(tubefaces), .true., 1, back=.true.)+1
                mergesepfaces = [tubefaces(ind2:), tubefaces(1:ind1)]
                sepfaces = tubefaces(ind1+1:ind2-1)
            else
                allocate(sepfaces(size(allhf)), mergesepfaces(size(alllf)))
                sepfaces = pack(tubefaces, ishfface(tubefaces))
                sepfaces = sepfaces(size(sepfaces):1:-1)
                mergesepfaces = pack(tubefaces, islfface(tubefaces))
            end if
        elseif (islfsep) then 
            if (ishfface(tubefaces(1)) .and. ishfface(tubefaces(size(tubefaces)))) then 
                ind1 = findloc(.not. ishfface(tubefaces), .true., 1, back=.false.)-1
                ind2 = findloc(.not. ishfface(tubefaces), .true., 1, back=.true.)+1
                mergesepfaces = [tubefaces(ind2:), tubefaces(1:ind1)]
                sepfaces = tubefaces(ind1+1:ind2-1)
            else
                allocate(sepfaces(size(alllf)), mergesepfaces(size(allhf)))
                sepfaces = pack(tubefaces, islfface(tubefaces))
                sepfaces = sepfaces(size(sepfaces):1:-1)
                mergesepfaces = pack(tubefaces, ishfface(tubefaces))
            end if
        end if         

        ! Extract sorted vertices from sorted merge separatrix faces
        call ExtractPolygonVertices(topomesh%face%vert(mergesepfaces, :), &
            size(mergesepfaces), mergevert)

        ! Sanity checks
        if (vert%type(mergevert(1)) == TMvertexsaddleID) then 
            print *, 'vertex: ', mergevert(1)
            print *, 'CollapseTMSeparatrixTubeTA: first sorted vertex of ' // & 
                'merging faces is X-point, unexpected. Aborting merge...'
            return 
        end if 
        if (vert%type(mergevert(size(mergevert))) == TMvertexsaddleID) then 
            print *, 'vertex: ', mergevert(size(mergevert))
            print *, 'CollapseTMSeparatrixTubeTA: last sorted vertex of ' // & 
                'merging faces is X-point, unexpected. Aborting merge...'
            return 
        end if 

        ! Insert X-points
        !================ 
        ! Sanity checks
        nxp = count(vert%type(mergevert) == TMvertexsaddleID)

        ! Sanity check
        if (nxp == 0) then
            print *, 'CollapseTMSeparatrixTubeTA: could not find any ' // & 
                'x-points, aborting merge...'
            return 
        end if

        ! Determine tracing direction
        if (ishfsep) then 
            ! Need to trace downstream, so opposite of the gradient
            tracedir = -1
        else
            ! Need to trace upstream, so along the gradient
            tracedir = 1
        end if 

        ! Determine lines
        allocate(xt(nxp), yt(nxp), xpointID(nxp), streamlines(nxp))
        xb(1) = -posinfval_R8()
        xb(2) = posinfval_R8()
        yb = xb
        cc = 0
        do i = 1, size(mergevert)
            if (vert%type(mergevert(i)) == TMvertexsaddleID) then
                ! Update counter
                cc = cc + 1

                ! Determine segment by finding the point that is closest
                ! to the x-point and of which the vector formed by the
                ! x-point and point lies within the domain (determined 
                ! by tangent vectors of separatrix faces neighbouring 
                ! the x-point)
                tf1 = mergesepfaces(i-1)
                tf2 = mergesepfaces(i)
                xpointID(cc) = mergevert(i)
                if (face%vert(tf1, 1) == mergevert(i)) then 
                    ! Need to take second vertex
                    x1 = face%x(tf1)%Get(2)
                    y1 = face%y(tf1)%Get(2)
                elseif  (face%vert(tf1, 2) == mergevert(i)) then 
                    ! Need to take N-1 vertex
                    x1 = face%x(tf1)%Get(face%x(tf1)%Size()-1)
                    y1 = face%y(tf1)%Get(face%y(tf1)%Size()-1)
                else
                    ! Unexpected
                    call gdErrorHandler('CollapseTMSeparatrixTubeTA: ' // &
                        'could not find x-point in face that is expected ' // & 
                        'to have x-point')
                end if
                if (face%vert(tf2, 1) == mergevert(i)) then 
                    ! Need to take second vertex
                    x2 = face%x(tf2)%Get(2)
                    y2 = face%y(tf2)%Get(2)
                elseif  (face%vert(tf2, 2) == mergevert(i)) then 
                    ! Need to take N-1 vertex
                    x2 = face%x(tf2)%Get(face%x(tf2)%Size()-1)
                    y2 = face%y(tf2)%Get(face%y(tf2)%Size()-1)
                else
                    ! Unexpected
                    print *, 'face vertices: ', face%vert(tf2, 1), face%vert(tf2, 2)
                    print *, 'face vertices: ', face%vert(tf1, 1), face%vert(tf1, 2)
                    call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                    call gdErrorHandler('CollapseTMSeparatrixTubeTA: ' // &
                        'could not find x-point in face that is expected ' // & 
                        'to have x-point')
                end if

                ! Construct the tangent vectors
                xxp = vert%x(mergevert(i))
                yxp = vert%y(mergevert(i))
                tx1 = x1 - xxp
                ty1 = y1 - yxp
                tx2 = x2 - xxp
                ty2 = y2 - yxp

                ! Determine the connecting point and face
                mindist = posinfval_R8()
                do j = 1, size(sepfaces)
                    ! Unpack
                    xs = face%x(sepfaces(j))%Get()
                    ys = face%y(sepfaces(j))%Get()

                    ! Compute (squared) distance
                    ds = (xs - xxp)**2 + (ys - yxp)**2

                    ! Compute vector
                    txv = xs - xxp
                    tyv = ys - yxp

                    ! Set distance to infinity where the vector does not
                    ! lie between the tangent vectors (i.e. this vector 
                    ! lies outside of the domain)
                    where ((txv*tx1 + tyv*ty1 < 0) .or. (txv*tx2 + tyv*ty2 < 0)) ds = posinfval_R8()

                    ! If all distances are infinite, skip
                    if (all(ds == posinfval_R8())) then 
                        cycle
                    end if 

                    ! Find the current minimal distance
                    tdist = minval(ds)

                    ! Check
                    if (tdist < mindist) then 
                        mindist = tdist
                        streamlines(cc)%x = [xxp, xs(minloc(ds, 1))]
                        streamlines(cc)%y = [yxp, ys(minloc(ds, 1))]
                    end if 

                end do
                
                ! Sanity check
                if (mindist == posinfval_R8()) then 
                    ! Not a single point found, unexpected
                    call gdErrorHandler('CollapseTMSeparatrixTubeTA: ' // & 
                        'could not find a single point that could ' // & 
                        'be connected to the x-point, unexpected ')
                end if 
            end if 
        end do 

        ! Determine streamline intersections
        nsl = size(streamlines)
        ntf = size(tubefaces)
        allocate(xintrda(nsl, ntf), yintrda(nsl, ntf), s1rrda(nsl, ntf), &
            s2rrda(nsl, ntf), s1ida(nsl, ntf), s2ida(nsl, ntf))
        do i = 1, nsl
            ! Unpack
            xs = streamlines(i)%x
            ys = streamlines(i)%y
            do j  = 1, ntf
                ! Compute intersections
                xf = face%x(tubefaces(j))%Get()
                yf = face%y(tubefaces(j))%Get()
                call SimplePolygonIntersections(xs, ys, xf, yf, xint, &
                    yint, s1, s2, s1r, s2r)

                ! Save intersections
                xintrda(i, j) = ConstructRealDynamicArray(xint)
                yintrda(i, j) = ConstructRealDynamicArray(yint)
                s1rrda(i, j) = ConstructRealDynamicArray(s1r)
                s2rrda(i, j) = ConstructRealDynamicArray(s2r)
                s1ida(i, j) = ConstructIntegerDynamicArray(s1)
                s2ida(i, j) = ConstructIntegerDynamicArray(s2)
                
            end do 
        end do 

        ! Process intersections & do checks
        allocate(newx(nsl), newy(nsl), sepfaceID(nsl), news2r(nsl), news2(nsl))
        do i = 1, nsl
            ! Check if one of the intersections is in a non-separatrix 
            ! boundary
            isintersectionfound = .false. 
            do j = 1, ntf
                if (any(tubefaces(j) == sepfaces)) then 
                    if (xintrda(i, j)%Size() /= 0) then 
                        ! Intersection found, need to do checks 
                        if (isintersectionfound) then 
                            ! Intersection was already found, so duplicate
                            ! intersections are present. Unexpected, abort merge
                            print *, 'CollapseTMSeparatrixTubeTA: multiple ' // &
                                'intersections found with separatrix, unexpected. Aborting merge'
                            return 
                        else
                            if (xintrda(i, j)%Size() > 1) then 
                                ! Intersection was already found, so duplicate
                                ! intersections are present. Unexpected, abort merge
                                print *, 'CollapseTMSeparatrixTubeTA: multiple ' // &
                                    'intersections found with separatrix with same face, unexpected. Aborting merge'
                                return 
                            else
                                ! Just one intersection, add that one
                                isintersectionfound = .true.
                                newx(i) = xintrda(i, j)%Get(1)
                                newy(i) = yintrda(i, j)%Get(1)
                                news2r(i) = s2rrda(i, j)%Get(1)
                                news2(i) = s2ida(i, j)%Get(1)
                                sepfaceID(i) = tubefaces(j)
                            end if 

                        end if 
                    end if 
                else
                    ! Can't have any intersections with other 
                    ! faces (yet) except merging faces
                    if (xintrda(i, j)%Size() /= 0) then 
                        if (any(tubefaces(j) == mergesepfaces)) then 
                            ! All good
                        else
                            print *, 'CollapseTMSeparatrixTubeTA: intersections ' // & 
                                'of streamline with non-separatrix face found, unexpected. Aborting merge...'
                            return 
                        end if 
                    end if 
                end if 
            end do 

            ! If we got here, check if an intersection was found
            if (.not. isintersectionfound) then 
                ! Exit
                print *, 'CollapseTMSeparatrixTubeTA: could not ' // &
                    'find intersection between streamline and separatrix, ' // & 
                    'exiting'
                return 
            end if 
        end do 

        ! Housekeeping
        end associate

        ! Determine flux surface ID for new faces
        call Unique(sepfaceID, usepfaceID)
        newfsID = topomesh%face%fsID(usepfaceID(1))
        do i = 1, size(usepfaceID)
            ! Sanity checks
            if (topomesh%face%fsID(usepfaceID(i)) /= newfsID) then 
                ! This can only happen if we didn't manage the flux surface
                ! IDs properly
                print *, 'CollapseTMSeparatixTubeTA: different flux surface ' // & 
                    'IDs found for same separatrix, tube construction in ' // & 
                    'grid generation stage may go wrong'
            end if 
        end do 

        ! Construct new aligned faces from old faces
        do i = 1, size(usepfaceID)
            ! Gather all intersections for this face
            nint = count(sepfaceID == usepfaceID(i))
            allocate(txp(nint), tnews2r(nint), tnews2(nint), sortindex(nint))
            txp = pack(xpointID, sepfaceID == usepfaceID(i))
            tnews2 = pack(news2, sepfaceID == usepfaceID(i))
            tnews2r = pack(news2r, sepfaceID == usepfaceID(i))

            ! Sort along s2r
            call Sort(tnews2r, ind=sortindex, ascend=.true.)
            txp = txp(sortindex)
            tnews2 = tnews2(sortindex)

            ! Add end points
            txp = [topomesh%face%vert(usepfaceID(i), 1), txp, topomesh%face%vert(usepfaceID(i), 2)]
            tnews2 = [1, tnews2, topomesh%face%pol(usepfaceID(i))%ne]

            ! New aligned faces are simply lines for now
            call ExtractTopologicalFacesFromPolygon(topomesh%face%pol(usepfaceID(i)), &
                txp, tnews2, topomesh%vert%x, topomesh%vert%y, v1, v2, xfnew, yfnew)

            ! Add to faces
            do j = 1, size(xfnew)
                if (topomesh%face%fsID(usepfaceID(i)) /= 0) then 
                    fsfval = topomesh%fsfval%Get(topomesh%face%fsID(usepfaceID(i)))
                else
                    fsfval = 0.0_R8
                end if 
                call AddTopologicalMeshFace(topomesh, [v1(j), v2(j)], xfnew(j), &
                    yfnew(j), topomesh%face%type(usepfaceID(i)), topomesh%face%fsID(usepfaceID(i)), fsfval)
            end do 

            ! Housekeeping
            deallocate(txp, tnews2r, tnews2, sortindex)
        end do 

        ! Update the adaptor data
        call tmadaptor%AddFaceData(topomesh)

        call WriteTopologicalMesh(topomesh, 'topomesh_temp1', .false.)

        ! Adjust topomesh
        !================
        ! Mark faces and vertices for removal
        allocate(remf(topomesh%face%ntot), remv(topomesh%vert%ntot))
        remf = .false. 
        remv = .false. 
        remf([mergesepfaces, usepfaceID]) = .true. ! separatrix faces to be removed
        where (topomesh%face%type(tuberadface) /= TMfacebndID) remf(tuberadface) = .true. ! remove non-boundary radial faces of tube
        do i = 1, size(mergevert)
            if (.not. any(topomesh%vert%type(mergevert(i)) == &  
                [TMvertexsaddleID, TMvertexbndID])) then 
                remv(mergevert(i)) = .true. ! boundary vertices will be removed by the simplification step
            end if 
        end do 

        ! Don't remove radial faces if the vertex also has a radial 
        ! face in the neighbouring tube (then also don't remove that vertex)
        do i = 1, size(mergevert)
            ! Get faces of vertex
            tvf = topomesh%vert%GetFace(mergevert(i))

            ! Check if it has a non-aligned face in the current tube 
            ! and a non-aligned face in the neighbouring tube. If so, 
            ! don't remove the vertex or it's non-aligned faces
            allocate(tvfna(0))
            do j = 1, size(tvf)
                if (any(topomesh%face%type(tvf(j)) == TMfacenonalignedID)) then 
                    tvfna = [tvfna, tvf(j)]
                end if 
            end do 
            if (any(istuberadface(tvfna)) .and. any(.not. istuberadface(tvfna))) then 
                ! Don't delete this vertex and faces
                remv(mergevert(i)) = .false.
                remf(tvfna) = .false.
            end if 

            ! Housekeeping
            deallocate(tvfna)
        end do

        ! Remove faces
        call RemoveTopologicalMeshFaceLogical(topomesh, remf)
        call tmadaptor%RemoveFaceDataLogical(remf)
        call WriteTopologicalMesh(topomesh, 'topomesh_temp2', .false.)

        ! Remove vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, remv)
        call WriteTopologicalMesh(topomesh, 'topomesh_temp3', .false.)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 
        call WriteTopologicalMesh(topomesh, 'topomesh_temp4', .false.)

        ! Split faces if necessary
        call SplitTopologicalMeshFaces(topomesh, remf)   
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 
        call WriteTopologicalMesh(topomesh, 'topomesh_temp5', .false.)
        ! Recompute all interconnections, cells, etc
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Data (recompute)
        call AddTopologicalMeshData(topomesh)

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Simplify separatrices
        !======================
        if (tmadaptor%allowsepsimplify) then 
            call tmadaptor%SimplifySeparatrices(topomesh)
        end if 

        ! Set output
        wasmerged = .true.
        
    end subroutine

    ! Separatrix simplification
    subroutine SimplifyTMSeparatricesTA(tmadaptor, topomesh)

        ! Description
        !============
        ! This routine simplifies the separatrix by applying a specific
        ! set of adaptations:
        ! - Merging of X-point(s) into one higher order X-point
        ! - Moving radial line points on X-points
        ! Note that these operations can in principle only occur after
        ! separatrix tubes have been merged or when connected topologies
        ! are encountered. Options should be propagated through the 
        ! topomesh adaptor object.

        ! It is *very* important to realize that this routine cannot
        ! hedge for any strange side effects that will occur when these
        ! adaptations are applied for X-points etc that lie far away 
        ! from each other, which is not the intended use. In fact, there 
        ! is likely no way to properly treat that, nor is it desired, 
        ! since these points are essential for the topological mesh and 
        ! since these should not pose any challenge during grid generation.  

        ! Algorithm
        !==========
        ! Simply check for each X-point each other vertex that is 
        ! connected through a separatrix segment. If the distance w.r.t.
        ! that vertex is lower than a defined tolerance, the vertex 
        ! it connects to is replaced by the X-point. This is done 
        ! iteratively until no more connections can be made. No checks
        ! are made a posteriori to see if faces cross etc, so handle 
        ! with care. The distance used here is the distance of the face.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor 
        type(TopomeshUDT), intent(inout)        :: topomesh

        ! Auxiliary
        integer(I8)                             :: tv 
        integer(I8), allocatable, dimension(:)  :: tf, tvf
        real(R8)                                :: lface
        real(R8), allocatable, dimension(:)     :: xf, yf
        logical, allocatable, dimension(:)      :: remf, remv
        type(RealDynamicArrayUDT)               :: xfrda, yfrda

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Loop until no more operations are done
        k = 1
        do while (k <= topomesh%vert%ntot)
            ! Cycle if this point is not an x-point
            if (topomesh%vert%type(k) /= TMvertexsaddleID) then 
                ! Increment counter
                k = k + 1

                ! Cycle
                cycle
            end if 

            ! Get all faces next to this X-point
            tf = topomesh%vert%GetFace(k)

            ! Loop over all faces
            do i = 1, size(tf)
                ! Check if it is a separatrix segment
                if (topomesh%face%type(tf(i)) /= TMfacesepID) then 
                    ! Skip
                    cycle
                end if
                
                ! Compute its distance
                xf = topomesh%face%x(tf(i))%Get()
                yf = topomesh%face%y(tf(i))%Get() 
                lface = sum(sqrt((xf(2:) - xf(1:size(xf)-1))**2 + &
                    (yf(2:) - yf(1:size(yf)-1))**2))

                ! Compare to minimal length
                if (lface >= tmadaptor%lpolsepmin) then 
                    ! Skip
                    cycle
                end if 

                ! Get non-x-point vertex
                if (topomesh%face%vert(tf(i), 1) == k) then 
                    tv = topomesh%face%vert(tf(i), 2)
                elseif (topomesh%face%vert(tf(i), 2) == k) then 
                    tv = topomesh%face%vert(tf(i), 1)
                else
                    ! Vertex not found, something wrong with topological
                    ! mesh interconnection
                    call gdErrorHandler('SimplifyTMSeparatrixTA: could ' // & 
                        'not find current X-point in face that should ' // & 
                        'have the X-point - unexpected error in topological ' // & 
                        'mesh interconnection')
                end if 

                ! Move the point: construct new faces and mark old ones
                ! for deletion
                allocate(remf(topomesh%face%ntot), remv(topomesh%vert%ntot))
                remf = .false. 
                remv = .false. 
                tvf = topomesh%vert%GetFace(tv)
                remf(tvf) = .true.
                remv(tv) = .true. 

                ! Construct and add new faces
                do j = 1, size(tvf)
                    ! Skip if the current face is the one being removed
                    if (tvf(j) == tf(i)) then 
                        cycle
                    end if 

                    ! Construct the new face
                    xfrda = topomesh%face%x(tvf(j))
                    yfrda = topomesh%face%y(tvf(j))
                    if (topomesh%face%vert(tvf(j), 1) == tv) then 
                        ! Adjust coordinates
                        call xfrda%Set(1, topomesh%vert%x(k))
                        call yfrda%Set(1, topomesh%vert%y(k))

                        ! Add face
                        call AddTopologicalMeshFace(topomesh, [k, topomesh%face%vert(tvf(j), 2)], xfrda, yfrda, &
                            topomesh%face%type(tvf(j)), topomesh%face%fsID(tvf(j)), 0.0_R8)
                    elseif (topomesh%face%vert(tvf(j), 2) == tv) then 
                        ! Adjust coordinates
                        call xfrda%Set(xfrda%Size(), topomesh%vert%x(k))
                        call yfrda%Set(yfrda%Size(), topomesh%vert%y(k))

                        ! Add face
                        call AddTopologicalMeshFace(topomesh, [topomesh%face%vert(tvf(j), 1), k], xfrda, yfrda, &
                            topomesh%face%type(tvf(j)), topomesh%face%fsID(tvf(j)), 0.0_R8)
                    else
                        call gdErrorHandler('SimplifyTMSeparatrixTA: ' // & 
                            'something wrong with interconnection')
                    end if 

                    ! Extend the deletion vector
                    remf = [remf, .false.]

                end do 

                ! Update the adaptor data
                call tmadaptor%AddFaceData(topomesh)

                ! Remove faces
                call RemoveTopologicalMeshFaceLogical(topomesh, remf)
                call tmadaptor%RemoveFaceDataLogical(remf)

                ! Remove vertices
                call RemoveTopologicalMeshVertexLogical(topomesh, remv)

                ! Simplify
                call SimplifyTopologicalMeshFaces(topomesh, remf)
                if (any(remf)) then
                    call tmadaptor%RemoveFaceData(remf)
                    call tmadaptor%AddFaceData(topomesh)
                end if 

                ! Split faces if necessary
                call SplitTopologicalMeshFaces(topomesh, remf)   
                if (any(remf)) then 
                    call tmadaptor%RemoveFaceData(remf)
                    call tmadaptor%AddFaceData(topomesh)
                end if 

                ! Recompute necessary interconnections
                ! Vertex faces
                call AddTopologicalMeshVertexFaces(topomesh)

                ! Reset k and cycle
                k = 1
                cycle
                
            end do 

            ! If we got here, then nothing happened. Update the counter
            k = k + 1
        end do 

        ! Update topomesh
        !================
        ! Data
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Data (recompute)
        call AddTopologicalMeshData(topomesh)

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

    end subroutine

    ! Tube splitting
    subroutine SplitTMTubesTA(tmadaptor, topomesh, tubes, dohfside, &
        dolfside, options)

        ! Description
        !============
        ! This routine splits tubes by inserting contours based on the
        ! merge criterion in order to form mergeable tubes for the 
        ! complex tube merging operator. It is checked whether the 
        ! tube should be split at high and/or low flux side. If possible, 
        ! at each desired side a contour is introduced at a distance so 
        ! that one of the merging criteria is almost exactly equal to the
        ! lower bound. This should result in a mergeable but non splittable
        ! tube. It is assumed this operation is only applied at tubes 
        ! that are a neighbour of tubes that do violate the merging 
        ! criterion. 

        ! Note: it is assumed that beforehand it is checked whether 
        ! tubes can be split. This is normally only the case if 
        ! the criterion value is at least twice the lower bound. This
        ! is checked here and an error is thrown if this is not 
        ! the case. 

        ! Modules
        !========
        use mod_search, only: findloc1D

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh
        integer(I8), dimension(:), intent(in)   :: tubes
        logical, dimension(:), intent(in)       :: dohfside, dolfside
        type(TopomeshOptionsUDT), intent(in)    :: options

        ! Auxiliary
        integer(I8)                             :: tf, ind
        integer(I8), allocatable, dimension(:)  :: tracetubes, tubeind, &
            tubef, temps1, temps2, tubefID, sortind, faceind, tracefaces, &
            tubecell, tempind
        real(R8)                                :: lffval, hffval, &
            psimin, psimax, lrad
        real(R8), allocatable, dimension(:)     :: hftracex, lftracex, &
            hftracey, lftracey, x, y, fval, dpsi, dlrad, dval, &
            thispsi, thislrad, dlc, tracex, tracey, temp, s2r, tempx, &
            tempy, temps1r, temps2r, xint, yint, newx, newy
        logical                                 :: isstartlf, wasfound
        logical, allocatable, dimension(:)      :: tracehf, tracelf, &
            tracec, keepind, remf, isincellpol
        type(ContourUDT), allocatable           :: allc(:), tempc(:)
        type(PolygonUDT)                        :: cellpol
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        type(RealDynamicArrayUDT)               :: xcrda, ycrda

        ! Loop
        integer(I8)                             :: i, j, k, cc

        ! Initialize
        !===========
        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_beforetubesplitting')
        end if 

        ! Initialize
        tracehf = dohfside 
        tracelf = dolfside 

        ! Check for trivial case
        if (.not. any(dohfside) .and. .not. any(dolfside)) then 
            return 
        end if 

        ! Compute 
        ! Note: include aligned boundaries here
        call tmadaptor%EvaluateTMTubesMergeCriterion(topomesh, includealbndin=.true.) 
        thispsi = tmadaptor%tubedpsi(tubes)
        thislrad = tmadaptor%tubelrad(tubes)
        dpsi = thispsi - tmadaptor%dpsimin 
        dlrad = thislrad - tmadaptor%lradmin

        ! Sanity checks
        dval = [thispsi - 3*tmadaptor%dpsimin, &
            thislrad - 3*tmadaptor%lradmin] ! if negative, then tube shouldn't have been marked for splitting
        if (any(dval < 0.0_R8)) then 
            call gdErrorHandler('SplitTMTubesTA: tubes ' // & 
                'were marked for splitting that are not wide enough, cannot continue')
        end if 

        ! Check for cases where both hf and lf should be traced if val > 3*lowerbound
        where (dohfside .and. dolfside .and. &
            ((thispsi < 4*tmadaptor%dpsimin) .or. (thislrad < 4*tmadaptor%lradmin))) tracelf = .false. 

        ! Set val$ues of criteria to inf where not active (i.e. dpsi = 0 
        ! or dlrad = 0)
        if (tmadaptor%dpsimin == 0.0_R8) then
            thispsi = posinfval_R8() 
            dpsi = posinfval_R8()
        end if 
        if (tmadaptor%lradmin == 0.0_R8) then
            thislrad = posinfval_R8() 
            dlrad = posinfval_R8()
        end if 

        ! Determine tracing points
        !=========================
        ! Initialize
        allocate(hftracex(size(tubes)), tracefaces(size(tubes)))
        hftracex = 0.0_R8
        lftracex = hftracex
        hftracey = hftracex
        lftracey = hftracex
        do i = 1, size(tubes)
            ! Get the tube radial faces
            tubef = topomesh%tube%GetFace(tubes(i))

            ! Check which criterion to follow and determine tracing points
            if (thispsi(i) < dlrad(i)) then ! Psi-based

                ! Here, we should be able to take any face of the tube, 
                ! as the bounds should be present on all radial faces
                tf = tubef(1)
                tracefaces(i) = tf

                ! Get the face coordinates and psi values
                x = topomesh%face%x(tf)%Get()
                y = topomesh%face%y(tf)%Get()
                fval = tmadaptor%facepsi_mono(tf)%Get()

                ! Check if psi increases or decreases
                isstartlf = (fval(size(fval)) - fval(1) >= 0.0_R8)

                ! Evaluate tube psi bounds - include aligned boundaries to prevent intersections with them
                call GetTMTubePsiLimitsTA(tmadaptor, topomesh, &
                    tubes(i), psimin, psimax, includealbndin=.true.)

                ! Compute flux values for tracing
                lffval = psimin + tmadaptor%dpsimin
                hffval = psimax - tmadaptor%dpsimin
                if (isstartlf) then 
                else
                    ! Switch fval, x, y for interpolation
                    fval = fval(size(fval):1:-1)
                    x = x(size(x):1:-1)
                    y = y(size(y):1:-1) 
                end if

                ! Determine tracing points by interpolation
                if (tracehf(i) .and. tracelf(i)) then 
                    ! Need to trace two contours
                    call Interpolate1D([lffval], temp, fval, x)
                    lftracex(i) = temp(1)
                    call Interpolate1D([lffval], temp, fval, y)
                    lftracey(i) = temp(1)
                    call Interpolate1D([hffval], temp, fval, x)
                    hftracex(i) = temp(1)
                    call Interpolate1D([hffval], temp, fval, y)
                    hftracey(i) = temp(1)
                    
                elseif (tracehf(i)) then 

                    call Interpolate1D([hffval], temp, fval, x)
                    hftracex(i) = temp(1)
                    call Interpolate1D([hffval], temp, fval, y)
                    hftracey(i) = temp(1)

                elseif (tracelf(i)) then 

                    call Interpolate1D([lffval], temp, fval, x)
                    lftracex(i) = temp(1)
                    call Interpolate1D([lffval], temp, fval, y)
                    lftracey(i) = temp(1)

                end if 

            else  ! Radial length based

                ! Re-evaluate the criterion and also query the face and
                ! the used length distribution
                call GetTMTubeRadialWidthTA(tmadaptor, topomesh, &
                    tubes(i), lrad, tf, dlc, includealbndin=.true.)
                tracefaces(i) = tf
                call GetTMTubePsiLimitsTA(tmadaptor, topomesh, &
                    tubes(i), psimin, psimax, includealbndin=.true.)

                ! Get the face coordinates and psi values
                x = topomesh%face%x(tf)%Get()
                y = topomesh%face%y(tf)%Get()
                fval = tmadaptor%facepsi_mono(tf)%Get()

                ! Check if psi increases or decreases
                isstartlf = (fval(size(fval)) - fval(1) >= 0.0_R8)

                ! Compute radial length values for tracing
                if (isstartlf) then 
                    lffval = tmadaptor%lradmin
                    hffval = dlc(size(x)) - tmadaptor%lradmin
                else
                    lffval = dlc(size(x)) - tmadaptor%lradmin
                    hffval = tmadaptor%lradmin
                end if

                ! Determine tracing points by interpolation
                if (tracehf(i) .and. tracelf(i)) then 
                    ! Need to trace two contours
                    call Interpolate1D([lffval], temp, dlc, x)
                    lftracex(i) = temp(1)
                    call Interpolate1D([lffval], temp, dlc, y)
                    lftracey(i) = temp(1)
                    call Interpolate1D([hffval], temp, dlc, x)
                    hftracex(i) = temp(1)
                    call Interpolate1D([hffval], temp, dlc, y)
                    hftracey(i) = temp(1)
                    
                elseif (tracehf(i)) then 

                    call Interpolate1D([hffval], temp, dlc, x)
                    hftracex(i) = temp(1)
                    call Interpolate1D([hffval], temp, dlc, y)
                    hftracey(i) = temp(1)

                elseif (tracelf(i)) then 

                    call Interpolate1D([lffval], temp, dlc, x)
                    lftracex(i) = temp(1)
                    call Interpolate1D([lffval], temp, dlc, y)
                    lftracey(i) = temp(1)

                end if 

            end if

            ! Check
            if (any(isnan([hftracex(i), hftracey(i), lftracex(i), lftracey(i)]))) then 
                print *,'SplitTMTubesTA: NaNs detected in tracing points'
            end if 
        end do 

        ! Trace contours
        !===============
        ! Concatenate for ease
        tracec = [tracehf, tracelf]
        tracex = [hftracex, lftracex]
        tracey = [hftracey, lftracey]
        tracetubes = [tubes, tubes]
        tracefaces = [tracefaces, tracefaces]

        ! Trace
        allocate(allc(0), tubeind(0), faceind(0))
        !$omp parallel do if (.not. omp_in_parallel()) &
        !$omp shared(tracec, tmadaptor, tracex, tracey, topomesh, allc, &
        !$omp tracetubes, tracefaces, tubeind, faceind) &
        !$omp private(tempc) default(none)
        do i = 1, size(tracec)
            if (tracec(i)) then 
                ! Trace the contour
                tempc = tmadaptor%fieldtracer%TraceContours([tracex(i)], [tracey(i)])
                
                ! Reformat into single contour
                if (size(tempc) == 1) then 
                    ! Do nothing, will add later on
                elseif (size(tempc) == 2) then 
                    ! Should be open contour
                    if (topomesh%tube%isclosed(tracetubes(i))) then 
                        call gdErrorHandler('SplitTMTubesTA: ' // & 
                            'tube is closed but contour is open, unexpected')
                    end if 

                    ! Concatenate 
                    tempc(1)%x = [tempc(2)%x(size(tempc(2)%x):2:-1), tempc(1)%x]
                    tempc(1)%y = [tempc(2)%y(size(tempc(2)%y):2:-1), tempc(1)%y]
                    tempc(1)%startsaddle = tempc(2)%endsaddle
                else
                    ! This may happen in very rare occasions if a 
                    ! separatrix segment was merged away and we are unlucky
                    ! enough to trace the contour at exactly the separatrix
                    ! value...
                    call gdErrorHandler('SplitTMTubesTA: ' // & 
                        'contour is either not found or consists of more than ' // & 
                        'two segments, unexpected')
                end if 

                ! Add
                !$omp critical
                allc = [allc, tempc(1)]
                tubeind = [tubeind, tracetubes(i)]
                faceind = [faceind, tracefaces(i)]
                !$omp end critical
            end if
        end do 
        !$omp end parallel do

        ! Process contours
        !=================
        ! Clean
        call CleanContours(allc)

        ! For open contours, check which parts to keep (only parts that 
        ! intersect with the tube faces)
        !$omp parallel do if (.not. omp_in_parallel()) & 
        !$omp default(none) &
        !$omp private(i, j, tubef, s2r, tubefID, xint, yint, tempx, &
        !$omp tempy, temps1, temps2, temps1r, temps2r, sortind, keepind, &
        !$omp ind, wasfound, tubecell, xcrda, ycrda, cellpol, cc, &
        !$omp isincellpol, newx, newy, tempind) &
        !$omp shared(topomesh, tubeind, allc, faceind)
        do i = 1, size(allc)
            ! Get tube faces
            tubef = topomesh%tube%GetFace(tubeind(i))

            ! Compute intersections
            allocate(s2r(0), tubefID(0), xint(0), yint(0))
            do j = 1, size(tubef)
                call SimplePolygonIntersections(topomesh%face%x(tubef(j))%Get(), &
                    topomesh%face%y(tubef(j))%Get(), allc(i)%x, allc(i)%y, &
                    tempx, tempy, temps1, temps2, temps1r, temps2r)

                ! Check if an intersection could not be found, or if multiple are
                ! found
                if (size(tempx) == 0) then 
                    ! Print warning - unexpected
                    print *, 'SplitTMTubesTA: traced ' // & 
                        'contour does not intersect with one of the radial ' // & 
                        'faces of tube ', tubeind(i), ', results may be unexpected'
                    print *, 'face vertices: ', topomesh%face%vert(tubef(j), :)
                    call Write2DCoordinateData(allc(i)%x, allc(i)%y, 'splitTMTubes_contour')
                elseif (size(tempx) > 1) then 
                    ! Print warning - unexpected
                    print *, 'SplitTMTubesTA: traced ' // & 
                        'contour intersects multiple times with one of the radial ' // & 
                        'faces of tube ', tubeind(i), ', results may be unexpected'
                end if

                ! Add intersection length coordinate to s2r
                s2r = [s2r, temps2r]
                xint = [xint, tempx]
                yint = [yint, tempy]
                tubefID = [tubefID, spread(tubef(j), 1, size(temps2r))]
            end do 

            ! Sort intersections
            allocate(sortind(size(s2r)))
            call Sort(s2r, ind=sortind, ascend=.true.)
            tubefID = tubefID(sortind)
            xint = xint(sortind)
            yint = yint(sortind)
            deallocate(sortind)

            ! Keep only parts that intersect with the tracing surface
            allocate(keepind(size(s2r)))
            keepind = .false. 
            if (allc(i)%isclosed .and. .not. topomesh%tube%isclosed(tubeind(i))) then
                ! closed contour for open tube - need to check differently. 
                ! Normally, the first and last intersection should be exactly
                ! in a tube face, since we trace from there. Since we
                ! only computed intersections with the current tube,
                ! there should be two segments that start at one tube 
                ! face and end in the other (errors are thrown otherwise)
                ! Of these two pieces, at least one should not lie 
                ! inside the tube polygon. This is checked by forming for
                ! one cell of the tube the cell polygon and checking if
                ! any vertices of the contour segment (except start/end) 
                ! are inside the cell polygon. 

                ! Checks
                if ((s2r(1) /= 0.0_R8) .or. (s2r(size(s2r)) /= size(allc(i)%x)-1)) then 
                    ! Unexpected
                    call gdErrorHandler('SplitTMTubesTA: '  // &
                        'contour does not seem to start and end in a tube face, ' // &
                        'unexpected since tracing should start from face')
                end if 

                ! Construct the cell polygon
                tubecell = topomesh%tube%GetCell(tubeind(i))
                call ConstructTopologicalMeshCellPolygon(topomesh, &
                    tubecell(1), xcrda, ycrda)
                call cellpol%Construct(xcrda%Get(), ycrda%Get())

                ! Loop over all segments
                cc = 0
                wasfound = .false. 
                do j = 1, size(s2r)-1
                    ! Check if this segment should be considered
                    if (((tubefID(j) == tubef(1)) .and. (tubefID(j+1) == tubef(2))) .or. &
                        ((tubefID(j) == tubef(2)) .and. (tubefID(j+1) == tubef(1)))) then 
                            
                        ! Update the counter
                        cc = cc + 1

                        ! Extract the segment coordinates without end points
                        tempind = [(k, k = 1, size(allc(i)%x))]
                        tempind = pack(tempind, (tempind > s2r(j)+1.0_R8) .and. &
                            (tempind < s2r(j+1)+1.0_R8))
                        tempx = [allc(i)%x(tempind)]
                        tempy = [allc(i)%y(tempind)]

                        ! Check sizes
                        if (size(tempx) == 0) then 
                            ! Need to refine a bit
                            print *, 'SplitTMTubesTA: closed contour segment ' // & 
                                'does not have any vertices in between, attempting ' // & 
                                'to refine (results may be inaccurate...)'
                            tempx = real([(k, k = 0, 3)], kind=R8)/3.0_R8*(xint(j+1) - xint(j)) + xint(j)
                            tempy = real([(k, k = 0, 3)], kind=R8)/3.0_R8*(yint(j+1) - yint(j)) + yint(j)
                            tempx = tempx(2:size(tempx)-1)
                            tempy = tempy(2:size(tempy)-1)
                        end if 

                        ! Check if there are any vertices in the cell polygon
                        call Inpolygon(cellpol, tempx, tempy, isincellpol)
                        if (any(isincellpol)) then 
                            ! Check if already found a segment
                            if (wasfound) then 
                                call gdErrorHandler('SplitTMTubesTA: ' // & 
                                    'both segments appear to be located ' // &
                                    'inside of tube, unexpected')
                            end if 
                            wasfound = .true.
                            newx = [xint(j), tempx, xint(j+1)] 
                            newy = [yint(j), tempy, yint(j+1)]
                        end if 
                    end if
                end do 

                ! Sanity check: should've found two segments
                if (cc /= 2) then 
                    call gdErrorHandler('SplitTMTubesTA: closed contour ' // & 
                        'did not have two segments as expected, this may be a bug')
                end if 

                ! Sanity check: should've found a segment that lies within
                ! the polygon cell
                if (.not. wasfound) then 
                    call gdErrorHandler('SplitTMTubesTA: could not find ' // & 
                        'a closed contour segment that lies within the tube')
                end if 


                ! Keep only part inbetween intersections
                allc(i)%x = newx
                allc(i)%y = newy

            elseif (allc(i)%isclosed .and. topomesh%tube%isclosed(tubeind(i))) then 
                ! Should be fine, nothing to check

            else ! both are open
                ! Check which ones to keep
                ind = findloc1D(tubefID, tubef)
                if (ind == 0) then 
                    ! Check for reverse order
                    ind = findloc1D(tubefID, tubef(size(tubef):1:-1))
                end if 
                if (ind == 0) then 
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call Write2DCoordinateData(allc(i)%x, allc(i)%y, 'splitTMTubes_contour')
                    call gdErrorHandler('SplitTMTubesTA: ' // & 
                        'could not find contour part that intersects with ' // & 
                        'all tube faces in the right order')
                end if 
                keepind(ind:ind+size(tubef)-1) = .true.

                ! Remove
                s2r = pack(s2r, keepind)
                xint = pack(xint, keepind)
                yint = pack(yint, keepind)

                ! Keep only part inbetween intersections
                allc(i)%x = [xint(minloc(s2r, 1)), &
                    allc(i)%x(ceiling(minval(s2r))+1:floor(maxval(s2r))+1), &
                    xint(maxloc(s2r, 1))]
                allc(i)%y = [yint(minloc(s2r, 1)), &
                    allc(i)%y(ceiling(minval(s2r))+1:floor(maxval(s2r))+1), &
                    yint(maxloc(s2r, 1))]

            end if 
            

            ! Housekeeping
            deallocate(s2r, xint, yint, tubefID, keepind)
        end do 
        !$omp end parallel do

        ! Clean again
        call CleanContours(allc)

        ! Add contours 
        !=============
        if (options%writedebugoutput) then 
            allocate(allpc(size(allc)))
            do i = 1, size(allc)
                call allpc(i)%Construct(allc(i)%x, allc(i)%y)
            end do 
            call tempps%Construct(allpc)
            call tempps%WriteData('splitcontours_all_beforeinsertion')
        end if 
        do i = 1, size(allc)
            ! Insert
            call InsertTopologicalMeshContour(topomesh, tmadaptor%magneticField, &
                allc(i), TMfacepolID, topomesh%nFs + i, remf)

            ! Update adaptor data
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)

            ! Set this flux surface as illegal surface to merge over 
            tmadaptor%splitfsIDs = [tmadaptor%splitfsIDs, topomesh%nFs + i]
        end do 
        topomesh%nFs = topomesh%nFs + size(allc)

        ! Rebuild topomesh
        !=================
        ! Trim the topological mesh
        call TrimTopologicalMesh(topomesh, tmadaptor%magneticField, &
            tmadaptor%vessel, remf)
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh, remf)
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Split
        call SplitTopologicalMeshFaces(topomesh, remf)  
        if (any(remf)) then 
            call tmadaptor%RemoveFaceData(remf)
            call tmadaptor%AddFaceData(topomesh)
        end if 

        ! Add necessary data
        !===================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)        

        ! Compute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_aftertubesplitting')
        end if

    end subroutine

    ! Clean-up of split TM tubes
    subroutine CleanSplitTMTubesTA(tmadaptor, topomesh)

        ! Description
        !============
        ! This routine serves to clean up any remaining splitted tubes
        ! that have not been involved in tube merging after the entire
        ! merging operation is done. In practice, we simply set the 
        ! merging criterion to infinity and call only the simple merging
        ! routine. This should only allow simple merges to be performed
        ! without any additional topological changes. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(inout)        :: topomesh 

        ! Auxiliary
        real(R8)                                :: lradmin_orig 
        logical                                 :: wasmerged

        ! Initialize
        !===========
        ! Copy original option for merging
        lradmin_orig = tmadaptor%lradmin

        ! Set length to infinity
        tmadaptor%lradmin = posinfval_R8()

        ! Simplify
        !=========
        ! Keep merging until no merge could be performed anymore
        do while (.true.)
            call tmadaptor%MergeTMTubesSimple(topomesh, wasmerged)
            if (.not. wasmerged) then 
                exit
            end if 
        end do 

        ! Reset
        !======
        tmadaptor%lradmin = lradmin_orig

    end subroutine

    ! Legal tube pair checking for merging
    function IsTubePairMergeableTA(tmadaptor, topomesh, hftube, lftube) &
        result(ismergeable)

        ! Description
        !============
        ! Check if the tubes on high flux and low flux side are 
        ! mergeable. A pair is considered to be mergeable if none of the
        ! aligned flux surfaces that they have in common have flux surface
        ! IDs that are prohibited to merge over. These are stored in the 
        ! tmadaptor object. Furthermore, the merge should not result in 
        ! a tube with overlapping psi values. Finally, any user-specified
        ! flux surface types that are prohibited (e.g. separatrix etc)
        ! are checked. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(in)           :: topomesh 
        integer(I8), dimension(:), intent(in)   :: hftube, lftube 
        logical                                 :: ismergeable 

        ! Auxiliary
        real(R8)                                :: newpsimax, newpsimin, &
            psimax, psimin
        integer(I8), allocatable, dimension(:)  :: mergefaces, &
            coreIDs, tc
        logical, allocatable, dimension(:)      :: isillegalface, &
            isillegalfsID, iscorecell

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Initialize
        ismergeable = .false. ! Will be set to true if all checks pass

        ! Unpack
        associate(&
            illegalfsIDs    => tmadaptor%illegalfsIDs,  &
            splitfsIDs      => tmadaptor%splitfsIDs     &
            )

        ! Check overlapping psi values
        !=============================
        ! Initialize
        allocate(isillegalface(topomesh%face%ntot), isillegalfsID(topomesh%nfs))
        isillegalfsID = .false.
        isillegalfsID(illegalfsIDs) = .true.
        isillegalfsID(splitfsIDs)   = .true.

        ! Set illegal faces to false if they don't have a flux surface ID
        where (topomesh%face%fsID /= 0)
            isillegalface = isillegalfsID(topomesh%face%fsID)
        elsewhere
            isillegalface = .false. 
        end where

        ! Check psi values
        newpsimax = posinfval_R8()
        newpsimin = -posinfval_R8()
        do i = 1, size(hftube)
            call GetTMTubePsiLimitsTA(tmadaptor, topomesh, hftube(i), psimin, psimax)
            newpsimax = min(newpsimax, psimax)
            if (psimin >= psimax) then 
                ! Not mergeable - return
                return
            end if 
        end do 
        do i = 1, size(lftube)            
            call GetTMTubePsiLimitsTA(tmadaptor, topomesh, lftube(i), psimin, psimax)
            newpsimin = max(newpsimin, psimin)
            if (psimin >= psimax) then 
                ! Not mergeable - return
                return
            end if 
        end do 
        if ((newpsimax - newpsimin) <= 0.0_R8) then 
            ! Not mergeable - return
            return
        end if 

        ! Check illegal flux surfaces
        !============================
        allocate(mergefaces(0))
        do i = 1, size(hftube)
            mergefaces = [mergefaces, topomesh%tube%GetLowFluxBndFace(hftube(i))]
        end do 
        do i = 1, size(lftube)
            mergefaces = [mergefaces, topomesh%tube%GetHighFluxBndFace(lftube(i))]
        end do 
        if (any(isillegalface(mergefaces))) then 
            ! Not mergeable - return
            return
        end if 

        ! Check if we merge over a flux surface that is prohibited by
        ! the user (e.g. separatrix, core, PF, ...)
        if (.not. tmadaptor%allowsepmerge) then 
            if (any(topomesh%face%type(mergefaces) == TMfacesepID)) then 
                return 
            end if 
        end if
        if (.not. tmadaptor%allowcoremerge) then 
            ! Core surfaces are only defined when core regions are already 
            ! removed - this may not yet be the case. Therefore we 
            ! determine core regions here and don't allow merging if one
            ! of the tubes has a core. 

            ! Check for faces
            if (any(topomesh%face%type(mergefaces) == TMfacecoreID)) then 
                return 
            end if 

            ! Check for core regions
            coreIDs = topomesh%GetCoreCellIDs()
            allocate(iscorecell(topomesh%cell%ntot))
            iscorecell = .false.
            iscorecell(coreIDs) = .true. 

            ! Check first set of tubes
            do i  = 1, size(hftube)
                tc =  topomesh%tube%GetCell(hftube(i))
                if (any(iscorecell(tc))) then 
                    return 
                end if 
            end do 

            ! Check second set of tubes
            do i  = 1, size(lftube)
                tc =  topomesh%tube%GetCell(lftube(i))
                if (any(iscorecell(tc))) then 
                    return 
                end if 
            end do 
        end if 
        if (.not. tmadaptor%allowpfmerge) then 
            if (any(topomesh%face%type(mergefaces) == TMfacePFID)) then 
                return 
            end if 
        end if 

        ! If we got here, all checks were passed so we can merge
        ismergeable = .true.

        ! Housekeeping
        end associate   

    end function

    ! Legal tube checking for splitting
    function IsTubeSplittableTA(tmadaptor, topomesh) result(issplittable)

        ! Description
        !============
        ! This routine checks for each tube in the topomesh whether 
        ! splitting is allowed. Normally, this is the case when the 
        ! tube is wide enough and when therefore the merging criterium
        ! itself is not fulfilled. Additionally, it is checked whether
        ! the tube does not have any overlapping psi values (this time
        ! including any aligned boundary faces) - in that case, the 
        ! contour would always intersect with the boundaries of the tube
        ! and hence cannot be properly split. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)                   :: tmadaptor
        type(TopomeshUDT), intent(in)               :: topomesh
        logical, allocatable, dimension(:)          :: issplittable

        ! Auxiliary
        real(R8), allocatable, dimension(:)         :: dpsi, dlrad

        ! Check based on criterion
        !=========================
        ! Compute criteria
        call tmadaptor%EvaluateTMTubesMergeCriterion(topomesh, includealbndin=.true.) ! do include here to prevent splitting of too small tubes
        dpsi = tmadaptor%tubedpsi - tmadaptor%dpsimin 
        dlrad = tmadaptor%tubelrad - tmadaptor%lradmin

        ! Check
        issplittable = (dpsi > 2*tmadaptor%dpsimin) .and. (dlrad > 2*tmadaptor%lradmin)

    end function

    ! Updaters for adaptor
    subroutine RemoveFaceDataLogicalTA(tmadaptor, remf)

        ! Description
        !============
        ! Remove any face data from faces that are/were removed in the
        ! topomesh. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        logical, dimension(:), intent(in)       :: remf 

        ! Initialize
        !===========
        ! Checks
        if (size(remf) /= size(tmadaptor%facepsi)) then 
            call gdErrorHandler('RemoveFaceDataLogicalTA: size of remf ' // &
                'does not correspond to size of facepsi')
        end if 
        if (size(remf) /= size(tmadaptor%facepsi_mono)) then 
            call gdErrorHandler('RemoveFaceDataLogicalTA: size of remf ' // &
                'does not correspond to size of facepsi_mono')
        end if 
        if (size(remf) /= size(tmadaptor%facedlcrad)) then 
            call gdErrorHandler('RemoveFaceDataLogicalTA: size of remf ' // &
                'does not correspond to size of facepsi')
        end if 

        ! Remove
        !=======
        tmadaptor%facepsi = pack(tmadaptor%facepsi, .not. remf)
        tmadaptor%facepsi_mono = pack(tmadaptor%facepsi_mono, .not. remf)
        tmadaptor%facedlcrad = pack(tmadaptor%facedlcrad, .not. remf)

    end subroutine

    subroutine AddFaceDataTA(tmadaptor, topomesh)

        ! Description
        !============
        ! This routine adds face of new faces in the topomesh. It is 
        ! assumed that faces are only added (so if faces are also
        ! deleted from the topomesh, it is assumed that RemoveFaceData
        ! is called first to remove these). New faces are assumed to 
        ! be appended to the topomesh. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(in)           :: topomesh

        ! Auxiliary
        integer(I8)                             :: nforig, nnewf, tfID
        real(R8), allocatable, dimension(:)     :: temp, temp_orig
        type(RealDynamicArrayUDT), allocatable, dimension(:)    :: &
            newfacepsi, newfacedlcrad, newfacepsi_mono

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Initialize
        nforig = size(tmadaptor%facepsi) 
        nnewf = topomesh%face%ntot - nforig

        ! Check
        if (nnewf < 0) then 
            call gdErrorHandler('AddFaceDataTA: less topomesh faces than ' // & 
                'face data in adaptor - first call face removal routine ' // & 
                'before adding face data')
        end if 

        ! Initialize further
        allocate(newfacepsi(nnewf), newfacedlcrad(nnewf), newfacepsi_mono(nnewf))

        ! Add data
        !=========
        ! Face psi values and radial length
        !$omp parallel do default(none) if ((.not. omp_in_parallel()) .and. &
        !$omp (nnewf >= 2*omp_get_num_threads())) &
        !$omp private(i, tfID, temp, temp_orig) &
        !$omp shared(topomesh, tmadaptor, nforig, nnewf, &
        !$omp newfacepsi, newfacedlcrad, newfacepsi_mono)
        do i = 1, nnewf
            ! Set current face ID
            tfID = nforig + i

            ! Compute
            call GetTMFacePsiValueDistribution(topomesh, tmadaptor%fieldtracer, tfID, &
                temp, temp_orig)
            newfacepsi_mono(i) = ConstructRealDynamicArray(temp)
            newfacepsi(i) = ConstructRealDynamicArray(temp_orig)
            temp = GetTMFaceRadialLengthDistribution(topomesh, &
                tmadaptor%fieldtracer, tmadaptor%magneticField, tfID)
            newfacedlcrad(i) = ConstructRealDynamicArray(temp)
        end do 
        !$omp end parallel do
        
        ! Append
        tmadaptor%facepsi = [tmadaptor%facepsi, newfacepsi]
        tmadaptor%facepsi_mono = [tmadaptor%facepsi_mono, newfacepsi_mono]
        tmadaptor%facedlcrad = [tmadaptor%facedlcrad, newfacedlcrad]


    end subroutine

    ! Aligned vessel part insertion
    subroutine InsertAlignedVesselParts(topomesh, fieldtracer, magneticField, &
        vessel, options)

        ! Description
        !============
        ! This routine inserts field aligned vessel parts into the 
        ! topological mesh *after* the full topological mesh has been
        ! constructed. The reason we do this afterwards, is because it
        ! is simply too difficult to hedge for some cases beforehand 
        ! without knowing the topology. Also, it is not generally possible
        ! to construct a conforming topological mesh for grid generation
        ! with a general allowance of field aligned walls (here, field 
        ! aligned does not necessarily mean true alignment, but rather a 
        ! vessel part that is indicated to be aligned and hence should be 
        ! treated as a flux surface boundary - a simple example of a full 
        ! aligned wall immediately shows that there is not always a 
        ! conforming topological mesh to be constructed). Therefore, we 
        ! limit the use of insertion of aligned vessel parts to some 
        ! specific cases. Currently, the following cases are supported:
        ! 
        ! - type 2 tangency points: insertion of aligned parts near
        !   tangency point only
        !
        ! If at some point other cases come up, they may be added here.

        ! It is VERY IMPORTANT to note that the introduction of 
        ! aligned vessel parts will break some assumptions/properties 
        ! of flux tubes:
        ! - tubes do not have unique high and low psi bounds anymore, 
        !   but rather a range of both. It still holds that the minimal
        !   high psi bound is always higher than the maximal low psi 
        !   bound. Therefore, this should have minor implications for
        !   grid generation etc
        ! - flux surfaces may intersect now with aligned vessel parts, 
        !   which may need to be hedged for in subsequent routines
        ! To adequately resolve these issues, we make use of the 
        ! aligned vessel part identifier for faces as defined in 
        ! mod_definitions. Vertices of aligned parts are retyped to 
        ! normal boundary vertices (TMvertexbndID), including tangency
        ! points. A new flux surface ID is made for the aligned 
        ! boundary and the flux surface value is set equal to some value
        ! that lies on that vessel part (this should be chosen properly 
        ! by the user), which should guarantee the above statements.
        
        ! Note 1: we're now hedging for deleted tangency points by checking
        ! when the sign of the angle between magnetic field and vessel
        ! changes. This indicates the presence of a (deleted) tangency 
        ! point, which may complicate contour tracing. In this case, the
        ! approximate location of the tangency point is taken as initial
        ! bound for the aligned vessel part. However, introducing a 
        ! contour there would likely lead to (again) very shallow 
        ! angles. Therefore, we override the case determination and set 
        ! it to case two, where the contour is traced in the middle of 
        ! the initial aligned vessel part. 
        
        ! Algorithms
        !===========
        ! Type 2 tangency points
        !-----------------------
        ! 1)    For each type 2 tangency point, determine which faces to 
        !       (partly) align and the alignment extent. It can only 
        !       be extended starting from the tangency point up to a 
        !       certain end point. Alignment is checked by comparing the
        !       magnetic field angle in each face edge center to the 
        !       minimal field angle defined in options%avpminangle (in 
        !       absolute value terms, avpminangle given in degrees). 
        !       Other criteria can of course also be implemented. If any
        !       faces were present that have to be extended, mark the 
        !       tangency point. 
        ! 2)    For each tube that has a vessel face to be (partly) 
        !       aligned, check in which of the following cases we are
        !       to determine how to adjust the tube and which parts to 
        !       insert (tangency point flux surface will always be removed):
        !   2.1)    The tube only has one face that is partly aligned. 
        !           Split the tube by inserting a contour at the 
        !           aligned/non-aligned transition. Keep track of the 
        !           contour flux surface ID - this will be used later 
        !           to determine which vessel parts are flux aligned 
        !           after contour insertion
        !   2.2)    The tube has one or two (more should be impossible)
        !           faces that are fully aligned. In this case, we 
        !           can't align both parts as then both tangency point
        !           faces should be deleted, which would yield a cell 
        !           that has more than two radial faces. Therefore, we 
        !           can only insert a contour somewhere along the tube
        !           (here, we take it at the middle of the first face
        !           - more optimal choices are possible as well)
        !   2.3)    The tube has two faces that are only partly aligned. 
        !           If no overlap: insert two contours at transition. If 
        !           overlap, but not the full tube is covered, insert a 
        !           single contour at the highest psi difference.
        !           Otherwise, do as in 2.2.
        !           Again, keep track of the flux surface IDs for later.
        ! 4)    Modify the contours: only the part that initially goes
        !       inside the vessel should be retained (until first 
        !       intersection with the vessel). Other contour parts 
        !       should not be added. This is effectively taken care of 
        !       by processing the contours somewhat analogously to 
        !       tangency point contours (keeping only the first part 
        !       until intersection with another vessel boundary) and 
        !       simplifying the topomesh afterwards. Additionally, we 
        !       need to hedge for discretization artefacts, i.e. the fact
        !       that traced contours might intersect with the original 
        !       contours due to differences in tracing resolution.
        ! 5)    Insert the contours, trim the topological mesh (remove 
        !       contours that are outside of domain) and simplify the
        !       topological mesh. Redefine face types by checking if
        !       a boundary face has a marked tangency point and a vertex
        !       with a newly introduced flux surface ID 
        ! 6)    Reconstruct tubes including aligned faces as potential 
        !       radial faces
        ! 7)    If a tube has an aligned face as radial face, then mark 
        !       the actual flux surface boundary of the tangency point 
        !       that was marked for removal. After all flux tubes have 
        !       been processed, delete these faces (and potentially any
        !       dangling vertices) and simplify the mesh. 
        ! 8)    Reconstruct all interconnections etc.

        ! Modules
        !========
        use mod_search, only: findloc1D
        use mod_sort, only: setdiff

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh
        class(ContourTracerUDT), intent(in)     :: fieldtracer 
        type(magneticFieldUDT), intent(in)      :: magneticField 
        type(TopomeshOptionsUDT), intent(in)    :: options 
        type(VesselUDT), intent(inout)          :: vessel

        ! Auxiliary
        logical                                 :: markface, ishftp(1:2), &
            hasnewfsIDnb1, hasnewfsIDnb2, wasfound
        logical, allocatable, dimension(:)      :: vertexmark, facemark, &
            isedgealigned, isentirefacealigned, remf, keepc, isstartface, &
            isalphapos, overridetubecase, remv, keepf, isincellpol
        integer(I8)                             :: thisf, tubecase, nfs, &
            ntpc, nint, nstc, intersectind, insertloc, tfmarktraceind(1:2)
        integer(I8), allocatable, dimension(:)  :: tf, tfbnd, afstartind, &
            afendind, tfmark, facevert, tfnb1, tfnb2, newfsIDs,  &
            markedtpIDs, sortind, vindI, vindJ, tsc, tfaceind, &
            vertexmarkIDs, tfnbv1, tfnbv2, tubecase_override, tvf, &
            tvmark, tubeID, tubecell, tubef, tempind, remtubeID, &
            difftubeID
        real(R8)                                :: avpminangle, tdl
        real(R8), allocatable, dimension(:)     :: tx, ty, xf, yf, dx, &
            dy, dn, bxf, byf, bnf, alpha, tpsinb1, tpsinb2, tpsitp, &
            xout, yout, iout, jout, tscr, dl, dlsum, thisx, thisy, &
            cosalpha, sinalpha, alphasigned, fval, dfval, tpsinew, &
            xint, yint, tempx, tempy, newx, newy

        type(ContourUDT), allocatable           :: tempc(:), allc(:)
        type(IntegerDynamicArrayUDT)            :: fsIDs, curvetypes, &
            cface, mergefsID
        type(IntegerDynamicArrayUDT), allocatable, dimension(:)     :: &
            sc, sf, faceind, contourind
        type(RealDynamicArrayUDT)               :: xcrda, ycrda
        type(RealDynamicArrayUDT), allocatable, dimension(:)        :: &
            sfr, scr, xintrda, yintrda
        type(PolygonUDT)                        :: cellpol
        type(PolygonUDT), allocatable           :: allpc(:)
        type(PolygonSetUDT)                     :: tempps
        
        ! Loop
        integer(I8)                             :: i, j, k, cc

        ! Initialize
        !===========
        ! Unpack 
        associate(&
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            tube        => topomesh%tube    &
            )

        ! Allocate
        allocate(vertexmark(vert%ntot), facemark(face%ntot), &
            afstartind(face%ntot), afendind(face%ntot), &
            isentirefacealigned(face%ntot), facevert(face%ntot), &
            tubecase_override(face%ntot), overridetubecase(face%ntot))

        ! Initialize
        nfs = topomesh%nFs
        vertexmark  = .false. 
        facemark    = .false.
        facevert    = 0
        afstartind  = 0
        afendind    = 0
        isentirefacealigned = .false.
        overridetubecase = .false. 
        tubecase_override = 0
        avpminangle = options%avpminangle/180.0_R8*pi_R8
        fsIDs       = ConstructIntegerDynamicArray()
        curvetypes  = ConstructIntegerDynamicArray()
        cface       = ConstructIntegerDynamicArray()
        mergefsID   = ConstructIntegerDynamicArray()

        ! Refine boundary faces if desired
        if (options%avprefinevessel) then 
            do i = 1, face%ntot 
                if (any(face%type(i) == [TMfacebndID, TMfacealbndID])) then 
                    call face%pol(i)%Refine(options%avpmaxvesseldist, &
                        options%avpminreffac)
                end if 
            end do 
        end if 

        ! Determine faces
        !================
        ! Tangency point type 1 faces
        do i = 1, vert%ntot
            if (vert%type(i) == TMvertextp2ID) then 
                ! Get faces
                tf = GetTMVertFace(vert, i)

                ! Get boundary faces 
                allocate(tfbnd(count(face%type(tf) == TMfacebndID))) 
                tfbnd = pack(tf, face%type(tf) == TMfacebndID)

                ! Sanity checks
                if (size(tfbnd) > 2) then 
                    ! This should not be possible
                    print *, 'vertex: ', i
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // & 
                        'tangency point has more than two vessel boundaries, ' // & 
                        'unexpected.')
                end if 

                ! Check for each face boundary face if it has parts to 
                ! align
                do j  = 1, size(tfbnd)
                    ! Initialize
                    markface = .false. 

                    ! Unpack
                    thisf   = tfbnd(j)

                    ! Compute metrics
                    tx = face%x(thisf)%Get()
                    ty = face%y(thisf)%Get()
                    dx = tx(2:) - tx(1:size(tx)-1)
                    dy = ty(2:) - ty(1:size(ty)-1)
                    dn = sqrt(dx**2 + dy**2)
                    dx = dx/dn 
                    dy = dy/dn

                    xf = 0.5*(tx(2:) + tx(1:size(tx)-1))
                    yf = 0.5*(ty(2:) + ty(1:size(ty)-1))
                    allocate(bxf(size(xf)), byf(size(xf)))
                    call magneticField%interp%Evaluate(xf, yf, 0, 1, bxf)
                    call magneticField%interp%Evaluate(xf, yf, 1, 0, byf)
                    bxf = -bxf
                    bnf = sqrt(bxf**2 + byf**2)
                    bxf = bxf/bnf 
                    byf = byf/bnf

                    ! Compute (absolute) angle 
                    cosalpha = bxf*dx + byf*dy
                    sinalpha = bxf*dy - byf*dx
                    alphasigned = atan2(sinalpha, cosalpha)
                    fval = fieldtracer%Evaluate(tx, ty)
                    dfval = fval(2:) - fval(1:size(fval)-1)
                    isalphapos = dfval > 0.0_R8
                    alpha = abs(alphasigned)
                    isedgealigned = (alpha < avpminangle) .or. ((pi_R8 - alpha) < avpminangle)                     
                    

                    ! Check if we should mark
                    if (any(isedgealigned)) then 
                        ! Still need to check if it starts at the 
                        ! tangency point
                        if (face%vert(thisf, 1) == i) then 
                            if (isedgealigned(1)) then 
                                ! Mark face and vertex
                                vertexmark(i)   = .true.
                                facemark(thisf) = .true.
                                facevert(thisf) = i

                                ! Compute face coordinate bound based on signed alpha + edge alignment
                                afstartind(thisf)   = 1 
                                if (all(isalphapos .eqv. isalphapos(1)) .and. &
                                    all(isedgealigned)) then 
                                    afendind(thisf) = size(tx)
                                    isentirefacealigned(thisf) = .true.
                                elseif (all(isalphapos .eqv. isalphapos(1))) then 
                                    ! Compute face coordinate bound based on aligned edges
                                    afendind(thisf) = findloc(isedgealigned, &
                                        .false., 1, back=.false.)
                                else
                                    ! Hit a removed tangency point 
                                    afendind(thisf) = findloc(isalphapos, &
                                        .not. isalphapos(1), 1, back=.false.)

                                    
                                    if (.not. all(isedgealigned)) then 
                                        ! Consider only edges til first non-aligned edge
                                        isedgealigned(findloc(isedgealigned, .false., 1, back=.false.):) = .false.
                                        ! if (findloc(isedgealigned, .false., 1, back=.false.) > afendind(thisf)) then 
                                        if (.false.) then 
                                            ! Bounded by alpha -> override to case 2
                                            afendind(thisf) = findloc((isalphapos .eqv. .not. isalphapos(1)) .or. &
                                                (isedgealigned), .true., 1, back=.true.) ! go as far as possible
                                            overridetubecase(thisf)     = .true.
                                            tubecase_override(thisf)    = 2 
                                        else
                                            ! Not bounded by alpha -> don't override and recompute
                                            afendind(thisf) = findloc(isedgealigned, .false., 1, back=.false.)
                                        end if 
                                    else
                                        ! Bounded by alpha -> override to case 2
                                        afendind(thisf) = findloc((isalphapos .eqv. .not. isalphapos(1)) .or. &
                                                (isedgealigned), .true., 1, back=.true.) ! go as far as possible
                                        overridetubecase(thisf)     = .true.
                                        tubecase_override(thisf)    = 2
                                    end if 
                                end if 
                            end if
                        elseif (face%vert(thisf, 2) == i) then 
                            if (isedgealigned(size(dx))) then 
                                ! Mark face and vertex
                                vertexmark(i)   = .true.
                                facemark(thisf) = .true.
                                facevert(thisf) = i

                                ! Compute face coordinate bound based on signed alpha
                                afendind(thisf)   = size(tx)
                                if (all(isalphapos .eqv. isalphapos(1)) .and. &
                                    all(isedgealigned)) then 
                                    afstartind(thisf) = 1
                                    isentirefacealigned(thisf) = .true.
                                elseif (all(isalphapos .eqv. isalphapos(1))) then 
                                    ! Compute face coordinate bound based on aligned edges
                                    afstartind(thisf) = findloc(isedgealigned, &
                                        .false., 1, back=.true.) + 1
                                else
                                    ! Hit a removed tangency point 
                                    afstartind(thisf) = findloc(isalphapos, &
                                        .not. isalphapos(size(tx)-1), 1, back=.true.) + 1
                                    if (.not. all(isedgealigned)) then 
                                        ! Consider only edges til first non-aligned edge
                                        isedgealigned(:findloc(isedgealigned, .false., 1, back=.true.)) = .false.
                                        ! if (findloc(isedgealigned, .false., 1, back=.true.)+1 < afendind(thisf)) then
                                        if (.false.) then  
                                            ! Bounded by alpha -> override to case 2
                                            afendind(thisf) = findloc((isalphapos .eqv. .not. isalphapos(size(tx)-1)) .or. &
                                                (isedgealigned), .true., 1, back=.false.) ! go as far as possible
                                            overridetubecase(thisf)     = .true.
                                            tubecase_override(thisf)    = 2 
                                        else
                                            ! Not bounded by alpha -> don't override and recompute
                                            afstartind(thisf) = findloc(isedgealigned, .false., 1, back=.true.) + 1
                                        end if 
                                    else
                                        ! Bounded by alpha -> override to case 2
                                        afendind(thisf) = findloc((isalphapos .eqv. .not. isalphapos(size(tx)-1)) .or. &
                                                (isedgealigned), .true., 1, back=.false.) ! go as far as possible
                                        overridetubecase(thisf)     = .true.
                                        tubecase_override(thisf)    = 2
                                    end if 
                                end if 
                            end if 
                        else
                            ! This shouldn't happen
                            print *, 'vertex: ', i
                            call WriteTopologicalMesh(topomesh, 'topomesh_error')
                            call gdErrorHandler('InsertAlignedVesselParts: ' // & 
                                'face of vertex does not contain vertex, unexpected')
                        end if 
                    end if 

                    ! Housekeeping
                    deallocate(bxf, byf)
                end do

                ! Housekeeping
                deallocate(tfbnd)
            end if 
        end do

        ! Compute marked vertex IDs for later
        allocate(vertexmarkIDs(count(vertexmark)))
        vertexmarkIDs = pack([(i, i = 1, topomesh%vert%ntot)], vertexmark)

        ! Trace contours
        !===============
        ! Trace depending on tube case
        allocate(allc(0), tubeID(0))
        do i = 1, tube%ntot 
            ! Get tube radial faces
            tf = GetTMTubeFace(tube, i)

            ! Get tube boundary faces
            if (allocated(tfbnd)) deallocate(tfbnd)
            allocate(tfbnd(count(face%type(tf) == TMfacebndID)))
            tfbnd = pack(tf, face%type(tf) == TMfacebndID)

            ! Determine tube hf/lf boundary faces 
            tfnb1 = GetTMTubeBndFace(tube, i, 1)
            tfnb2 = GetTMTubeBndFace(tube, i, 2)

            ! Sanity checks
            if (size(tfbnd) > 2) then 
                ! This is not supported
                print *, 'tube: ', i
                call WriteTopologicalMesh(topomesh, 'topomesh_error')
                call gdErrorHandler('InsertAlignedVesselParts: tube ' // &
                    'has more than two boundary faces, unexpected. Cannot proceed.')
            elseif (size(tfbnd) == 0) then 
                ! No boundary faces -> nothing to do here (probably 
                ! closed tube in a core region or similar)
                cycle
            elseif (size(tfbnd) == 1) then 
                ! This is not supported
                print *, 'tube: ', i
                call WriteTopologicalMesh(topomesh, 'topomesh_error')
                call gdErrorHandler('InsertAlignedVesselParts: tube ' // &
                    'has only one boundary face, unexpected. Cannot proceed.')
            end if

            ! Determine which case we're in 
            tubecase = 0
            if (allocated(tfmark)) deallocate(tfmark)
            allocate(tfmark(count(facemark(tfbnd))))
            tfmark = pack(tfbnd, facemark(tfbnd))
            if (size(tfmark) == 0) then 
                ! Nothing to do here
                cycle 
            elseif (size(tfmark) == 1) then 
                ! Single face -> need to check bounds
                if (isentirefacealigned(tfmark(1))) then 
                    tubecase = 2
                else
                    tubecase = 1
                end if 
            else
                ! Two boundary faces - checked before
                if (any(isentirefacealigned(tfmark))) then 
                    tubecase = 2
                else
                    tubecase = 3
                end if 
            end if 

            ! Get vertices
            tvmark = [topomesh%face%vert(tfmark, 1), topomesh%face%vert(tfmark, 2)]
            tvmark = pack(tvmark, vertexmark(tvmark))

            ! Override 
            if (any(overridetubecase(tfmark))) then 
                ! Just take one...
                tubecase = tubecase_override(tfmark(1))
            end if 

            ! Trace and insert contours where necessary
            select case (tubecase)

            case (1) 
                
                ! One face, partly aligned -> insert contour
                if (face%vert(tfmark(1), 1) == facevert(tfmark(1))) then 
                    ! First vertex is start -> trace contour at transition
                    tempc = fieldtracer%TraceContours(&
                        [face%x(tfmark(1))%Get(afendind(tfmark(1)))], &
                        [face%y(tfmark(1))%Get(afendind(tfmark(1)))])
                    allc = [allc, tempc]
                    nfs = nfs + 1
                    call fsIDs%Append(spread(nfs, 1, size(tempc)))
                    call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                    call cface%Append(spread(tfmark(1), 1, size(tempc)))

                    ! Add tube ID of contours
                    tubeID = [tubeID, spread(i, 1, size(tempc))]
                else
                    ! Second vertex is start
                    tempc = fieldtracer%TraceContours(&
                        [face%x(tfmark(1))%Get(afstartind(tfmark(1)))], &
                        [face%y(tfmark(1))%Get(afstartind(tfmark(1)))])
                    allc = [allc, tempc]
                    nfs = nfs + 1
                    call fsIDs%Append(spread(nfs, 1, size(tempc)))
                    call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                    call cface%Append(spread(tfmark(1), 1, size(tempc)))

                    ! Add tube ID of contours
                    tubeID = [tubeID, spread(i, 1, size(tempc))]
                end if 

                

            case (2) 
                
                ! One or two faces, fully aligned -> insert contour
                ! close to middle of first face (needs to start in 
                ! face point...)

                ! Compute metrics
                tx = face%x(tfmark(1))%Get()
                ty = face%y(tfmark(1))%Get()
                dx = tx(2:) - tx(1:size(tx)-1)
                dy = ty(2:) - ty(1:size(ty)-1)
                dl = sqrt(dx**2 + dy**2)
                dlsum = tx
                dlsum(1) = 0.0_R8 
                do j = 2, size(dlsum)
                    dlsum(j) = dlsum(j-1) + dl(j-1)
                end do 
                
                ! Normalize
                dlsum = dlsum/dlsum(size(dlsum))

                ! Compute tracing point
                tdl = dlsum(afstartind(tfmark(1)))+0.5_R8*(dlsum(afendind(tfmark(1))) - dlsum(afstartind(tfmark(1))))
                call Interpolate1D([tdl], thisx, dlsum, tx)
                call Interpolate1D([tdl], thisy, dlsum, ty)

                ! Insert tracing point in face and reconstruct polygon
                insertloc = findloc(abs(dlsum - tdl)*sum(dl) < 1e-6, .true., 1, back=.false.)
                if ((insertloc /= 0) .and. (insertloc /= 1) .and. (insertloc /= size(tx))) then 
                    ! Relocate existing point
                    call face%x(tfmark(1))%Set([insertloc], thisx)
                    call face%y(tfmark(1))%Set([insertloc], thisy)
                else 
                    ! Insert in the middle
                    insertloc = findloc(dlsum > tdl, .true., 1, back=.false.)
                    if (insertloc == 0) then 
                        ! This should never happen, unless the original face length 
                        ! is near machine precision and round-off effects 
                        ! come in play
                        call gdErrorHandler('InsertAlignedVesselParts: could ' // & 
                            'not insert point in face, unexpected')
                    end if

                    ! Insert
                    call face%x(tfmark(1))%Insert(thisx, [insertloc])
                    call face%y(tfmark(1))%Insert(thisy, [insertloc])
                end if 
                    
                ! Reconstruct polygon
                call face%pol(tfmark(1))%Construct(face%x(tfmark(1))%Get(), &
                    face%y(tfmark(1))%Get())

                ! Trace contour and add
                tempc = fieldtracer%TraceContours(thisx, thisy)
                allc = [allc, tempc]
                nfs = nfs + 1
                call fsIDs%Append(spread(nfs, 1, size(tempc)))
                call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                call cface%Append(spread(tfmark(1), 1, size(tempc)))

                ! Add tube ID of contours
                tubeID = [tubeID, spread(i, 1, size(tempc))]


            case (3)

                ! Two faces that are partially aligned, need to check further

                ! Sanity check: should have two marked vertices
                if (size(tvmark) /= 2) then 
                    print *, 'found marked vertices: ', tvmark
                    call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.) 
                    call gdErrorHandler('InsertAlignedVesselParts: could ' // &
                        'not find exactly two marked vertices, unexpected')
                end if 

                ! Get tube and tangency points psi values
                tfnb1 = GetTMTubeHighFluxBndFace(topomesh%tube, i)
                tfnb2 = GetTMTubeLowFluxBndFace(topomesh%tube, i)
                tpsinb1 = topomesh%fsfval%Get(face%fsID(tfnb1))
                tpsinb2 = topomesh%fsfval%Get(face%fsID(tfnb2))
                tpsitp = topomesh%fsfval%Get(tvmark)

                ! Compute psi values of potential new contours and their
                ! tracing index
                if (vertexmark(face%vert(tfmark(1), 1))) then 
                    ! First vertex is tangency point, so need to look at 
                    ! end index
                    tfmarktraceind(1) = afendind(tfmark(1))
                    
                elseif (vertexmark(face%vert(tfmark(1), 2))) then 
                    ! Second vertex is tangency point, so need to look 
                    ! at start index
                    tfmarktraceind(1) = afstartind(tfmark(1))

                else
                    ! No vertices are tangency points - unexpected
                    call gdErrorHandler('InsertAlignedVesselParts: no ' // & 
                        'marked tangency points found in marked face')
                end if 
                if (vertexmark(face%vert(tfmark(2), 1))) then 
                    ! First vertex is tangency point, so need to look at 
                    ! end index
                    tfmarktraceind(2) = afendind(tfmark(2))
                    
                elseif (vertexmark(face%vert(tfmark(2), 2))) then 
                    ! Second vertex is tangency point, so need to look 
                    ! at start index
                    tfmarktraceind(2) = afstartind(tfmark(2))

                else
                    ! No vertices are tangency points - unexpected
                    call gdErrorHandler('InsertAlignedVesselParts: no ' // & 
                        'marked tangency points found in marked face')
                end if 
                tpsinew = [fieldtracer%Evaluate([face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                    [face%y(tfmark(1))%Get(tfmarktraceind(1))]), &
                fieldtracer%Evaluate([face%x(tfmark(2))%Get(tfmarktraceind(2))], &
                    [face%y(tfmark(2))%Get(tfmarktraceind(2))])]

                ! Determine if first tp is high/low field tp
                if (any([face%vert(tfnb1, 1), face%vert(tfnb1, 2)] == tvmark(1))) then 
                    ishftp(1) = .true. 
                elseif (any([face%vert(tfnb2, 1), face%vert(tfnb2, 2)] == tvmark(1))) then 
                    ishftp(1) = .false. 
                else
                    ! Shouldn't happen
                    print *, 'tube: ', i, 'vertex: ', facevert(tfmark(1))
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // &
                        'tangency point not found in any tube aligned boundary')
                end if 

                ! Determine if second tp is high/low field tp
                if (any([face%vert(tfnb1, 1), face%vert(tfnb1, 2)] == tvmark(2))) then 
                    ishftp(2) = .true. 
                elseif (any([face%vert(tfnb2, 1), face%vert(tfnb2, 2)] == tvmark(2))) then 
                    ishftp(2) = .false. 
                else
                    ! Shouldn't happen
                    print *, 'tube: ', i, 'vertex: ', facevert(tfmark(2))
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // &
                        'tangency point not found in any tube aligned boundary')
                end if 
                
                ! Check which contours to trace
                if (all(ishftp)) then 
                    ! Take lowest psi value
                    if (tpsitp(1) < tpsitp(2)) then 
                        tempc = fieldtracer%TraceContours(&
                            [face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                            [face%y(tfmark(1))%Get(tfmarktraceind(1))])
                        allc = [allc, tempc]
                        nfs = nfs + 1
                        call fsIDs%Append(spread(nfs, 1, size(tempc)))
                        call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                        call cface%Append(spread(tfmark(1), 1, size(tempc)))

                        ! Add tube ID of contours
                        tubeID = [tubeID, spread(i, 1, size(tempc))]
                    else
                        tempc = fieldtracer%TraceContours(&
                            [face%x(tfmark(2))%Get(tfmarktraceind(2))], &
                            [face%y(tfmark(2))%Get(tfmarktraceind(2))])
                        allc = [allc, tempc]
                        nfs = nfs + 1
                        call fsIDs%Append(spread(nfs, 1, size(tempc)))
                        call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                        call cface%Append(spread(tfmark(2), 1, size(tempc)))

                        ! Add tube ID of contours
                        tubeID = [tubeID, spread(i, 1, size(tempc))]
                    end if
                elseif (.not. any(ishftp)) then 
                    ! Take highest psi value
                    if (tpsitp(1) > tpsitp(2)) then 
                        tempc = fieldtracer%TraceContours(&
                            [face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                            [face%y(tfmark(1))%Get(tfmarktraceind(1))])
                        allc = [allc, tempc]
                        nfs = nfs + 1
                        call fsIDs%Append(spread(nfs, 1, size(tempc)))
                        call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                        call cface%Append(spread(tfmark(1), 1, size(tempc)))

                        ! Add tube ID of contours
                        tubeID = [tubeID, spread(i, 1, size(tempc))]
                    else
                        tempc = fieldtracer%TraceContours(&
                            [face%x(tfmark(2))%Get(tfmarktraceind(2))], &
                            [face%y(tfmark(2))%Get(tfmarktraceind(2))])
                        allc = [allc, tempc]
                        nfs = nfs + 1
                        call fsIDs%Append(spread(nfs, 1, size(tempc)))
                        call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                        call cface%Append(spread(tfmark(2), 1, size(tempc)))

                        ! Add tube ID of contours
                        tubeID = [tubeID, spread(i, 1, size(tempc))]
                    end if
                else
                    ! Lowest and highest, need to check overlap
                    if (ishftp(1)) then 
                        ! First is highest, second is lowest. 
                        if (tpsinew(1) < tpsinew(2)) then 
                            ! Overlap -> just take first one
                            !if (abs(tpsitp(1) - highpsi) < abs(tpsitp(2) - lowpsi)) then 
                                tempc = fieldtracer%TraceContours(&
                                    [face%x(tfmark(1))%Get(afendind(tfmark(1)))], &
                                    [face%y(tfmark(1))%Get(afendind(tfmark(1)))])
                                allc = [allc, tempc]
                                nfs = nfs + 1
                                call fsIDs%Append(spread(nfs, 1, size(tempc)))
                                call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                                call cface%Append(spread(tfmark(1), 1, size(tempc)))
                            !else
                            !    tempc = fieldtracer%TraceContours(&
                            !        [face%x(tfmark(2))%Get(afendind(tfmark(2)))], &
                            !        [face%y(tfmark(2))%Get(afendind(tfmark(2)))])
                            !    allc = [allc, tempc]
                            !    nfs = nfs + 1
                            !    call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            !    call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            !    call cface%Append(spread(tfmark(2), 1, size(tempc)))
                            !end if 
                        else 
                            ! No overlap -> trace both contours
                            tempc = fieldtracer%TraceContours(&
                                [face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                                [face%y(tfmark(1))%Get(tfmarktraceind(1))])
                            allc = [allc, tempc]
                            nfs = nfs + 1
                            call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            call cface%Append(spread(tfmark(1), 1, size(tempc)))

                            ! Add tube ID of contours
                            tubeID = [tubeID, spread(i, 1, size(tempc))]

                            tempc = fieldtracer%TraceContours(&
                                [face%x(tfmark(2))%Get(tfmarktraceind(2))], &
                                [face%y(tfmark(2))%Get(tfmarktraceind(2))])
                            allc = [allc, tempc]
                            nfs = nfs + 1
                            call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            call cface%Append(spread(tfmark(2), 1, size(tempc)))

                            ! Add tube ID of contours
                            tubeID = [tubeID, spread(i, 1, size(tempc))]
                        end if 
                    else
                        ! Second is highest
                        if (tpsitp(1) > tpsitp(2)) then 
                            ! Overlap -> take one with lowest delta Psi
                            !if (abs(tpsitp(1) - lowpsi) < abs(tpsitp(2) - highpsi)) then 
                                tempc = fieldtracer%TraceContours(&
                                    [face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                                    [face%y(tfmark(1))%Get(tfmarktraceind(1))])
                                allc = [allc, tempc]
                                nfs = nfs + 1
                                call fsIDs%Append(spread(nfs, 1, size(tempc)))
                                call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                                call cface%Append(spread(tfmark(1), 1, size(tempc)))
                            !else
                            !    tempc = fieldtracer%TraceContours(&
                            !        [face%x(tfmark(2))%Get(afendind(tfmark(2)))], &
                            !        [face%y(tfmark(2))%Get(afendind(tfmark(2)))])
                            !    allc = [allc, tempc]
                            !    nfs = nfs + 1
                            !    call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            !    call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            !    call cface%Append(spread(tfmark(2), 1, size(tempc)))
                            !end if 
                        else 
                            ! No overlap -> trace both contours
                            tempc = fieldtracer%TraceContours(&
                                [face%x(tfmark(1))%Get(tfmarktraceind(1))], &
                                [face%y(tfmark(1))%Get(tfmarktraceind(1))])
                            allc = [allc, tempc]
                            nfs = nfs + 1
                            call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            call cface%Append(spread(tfmark(1), 1, size(tempc)))

                            ! Add tube ID of contours
                            tubeID = [tubeID, spread(i, 1, size(tempc))]

                            tempc = fieldtracer%TraceContours(&
                                [face%x(tfmark(2))%Get(tfmarktraceind(2))], &
                                [face%y(tfmark(2))%Get(tfmarktraceind(2))])
                            allc = [allc, tempc]
                            nfs = nfs + 1
                            call fsIDs%Append(spread(nfs, 1, size(tempc)))
                            call curvetypes%Append(spread(TMfacepolID, 1, size(tempc)))
                            call cface%Append(spread(tfmark(2), 1, size(tempc)))

                            ! Add tube ID of contours
                            tubeID = [tubeID, spread(i, 1, size(tempc))]
                        end if 
                    end if 
                end if 

            end select

            ! Housekeeping
            deallocate(tfbnd, tfmark)
            
        end do

        ! Clean
        call CleanContours(allc)

        !! Retype fully aligned parts
        !do i = 1, face%ntot
        !    if (isentirefacealigned(i)) then 
        !        face%type(i) = TMfacealbndID 
        !    end if 
        !end do 

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp0')
        end if 

        ! Post-process
        !=============
        ! This part is quite related to tangency point contour processing
        ! as done in AddTopologicalMeshContours. We compute the 
        ! intersections with all boundary faces and trim the contour based 
        ! on those intersections. We do not adjust the original vessel 
        ! faces. 

        ! Initialize
        ntpc = size(allc)
        allocate(keepc(ntpc))
        keepc = .true. 
        allocate(sc(ntpc), scr(ntpc), faceind(ntpc), sf(topomesh%face%ntot), &
            sfr(topomesh%face%ntot), contourind(topomesh%face%ntot), &
            xintrda(ntpc), yintrda(ntpc))
        do i = 1, ntpc
            sc(i) = ConstructIntegerDynamicArray()
            scr(i) = ConstructRealDynamicArray()
            xintrda(i) = ConstructRealDynamicArray()
            yintrda(i) = ConstructRealDynamicArray()
            faceind(i) = ConstructIntegerDynamicArray()
        end do 
        do i = 1, topomesh%face%ntot
            sf(i) = ConstructIntegerDynamicArray()
            sfr(i) = ConstructRealDynamicArray()
            contourind(i) = ConstructIntegerDynamicArray()
        end do 

        ! Compute intersections
        !$omp parallel do default(none) &
        !$omp schedule(dynamic) & 
        !$omp private(xout, yout, vindI, vindJ, iout, jout, i, j, tubef) & 
        !$omp shared(ntpc, topomesh, allc, sc, scr, faceind, sf, sfr, &
        !$omp contourind, xintrda, yintrda, tubeID)
        do i = 1, ntpc
            ! Get all faces of this tube
            tubef = [topomesh%tube%GetBndFace(tubeID(i), 1), &
                topomesh%tube%GetBndFace(tubeID(i), 2), &
                topomesh%tube%GetFace(tubeID(i))]
            do k = 1, size(tubef)
                j = tubef(k)
                ! Compute intersections
                call SimplePolygonIntersections(allc(i)%x, allc(i)%y, &
                    topomesh%face%x(j)%Get(), topomesh%face%y(j)%Get(), &
                    xout, yout, vindI, vindJ, iout, jout)

                ! Add, if any
                if (allocated(xout)) then 
                    if (size(xout) > 0) then 
                        !$omp critical
                        ! Add intersection data to contour
                        call sc(i)%Append(vindI)
                        call scr(i)%Append(iout)
                        call xintrda(i)%Append(xout)
                        call yintrda(i)%Append(yout)
                        call faceind(i)%Append(spread(j, 1, size(vindI)))

                        ! Add intersection data to face
                        call sf(j)%Append(vindJ)
                        call sfr(j)%Append(jout)
                        call contourind(j)%Append(spread(i, 1, size(vindI))) 
                        !$omp end critical
                    end if 
                end if
            end do 
        end do
        !$omp end parallel do 

        ! Process contours
        do i = 1, ntpc
            ! Unpack
            tsc = sc(i)%Get()
            tfaceind = faceind(i)%Get()
            tscr = scr(i)%Get()
            xint = xintrda(i)%Get()
            yint = yintrda(i)%Get()
            nint = size(tsc)
            nstc = size(allc(i)%x)-1

            ! Sort 
            sortind = [(k, k = 1, nint)]
            call Sort(tscr, ind=sortind, ascend=.true.)
            tsc = tsc(sortind)
            xint = xint(sortind)
            yint = yint(sortind)
            tfaceind = tfaceind(sortind)

            ! Check
            if (tscr(1) /= 0.0_R8) then 
                ! First intersection should always be in tangency point
                call gdErrorHandler('InsertAlignedVesselParts: first ' // & 
                    'intersection of contour is not in boundary' )
            end if 
            if (allc(i)%isclosed) then 
                ! Last intersection should also be in tangency point
                if (tscr(nint) /= size(allc(i)%x)-1) then 
                    call gdErrorHandler('InsertAlignedVesselParts: last ' // & 
                    'intersection of closed contour is not in boundary' )
                end if 
            end if 

            ! Check if there are any intersections with an aligned boundary
            ! -> should not add contour as it intersects tube boundaries
            if (size(GetCommonElements(topomesh%face%type(tfaceind), TMfacealignedID)) /= 0) then 
                keepc(i) = .false. 
                cycle
            end if

            ! Determine which intersections were in the starting face
            isstartface = tfaceind == cface%Get(i)   

            ! Get the tube start and end face
            tf = topomesh%tube%GetFace(tubeID(i))
            tubef = [tf(1), tf(size(tf))]

            ! Determine which part(s) of contour(s) to keep
            if (allc(i)%isclosed) then ! closed contour
                ! closed contour for open tube - need to check differently. 
                ! Normally, the first and last intersection should be exactly
                ! in a tube face, since we trace from there. Since we
                ! only computed intersections with the current tube,
                ! there should be two segments that start at one tube 
                ! face and end in the other (errors are thrown otherwise)
                ! Of these two pieces, at least one should not lie 
                ! inside the tube polygon. This is checked by forming for
                ! one cell of the tube the cell polygon and checking if
                ! any vertices of the contour segment (except start/end) 
                ! are inside the cell polygon. 
                
                ! Print message
                print *, 'InsertAlignedVesselParts: closed contour for ' // & 
                    'open tube detected, code not yet verified'

                ! Keep only intersections that intersect in one of the 
                ! tube's bounding radial faces
                allocate(keepf(size(tfaceind)))
                keepf = .false. 
                where (tfaceind == tubef(1) .or. tfaceind == tubef(2)) keepf = .true. 
                tscr = pack(tscr, keepf)
                xint = pack(xint, keepf)
                yint = pack(yint, keepf)
                tfaceind = pack(tfaceind, keepf)
                deallocate(keepf)

                ! Sanity check
                if (all(isstartface)) then 
                    ! No intersections with other boundaries - this should
                    ! not occur!
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // & 
                        'detected closed contour that only intersects in ' // & 
                        'starting face - unexpected')
                end if 

                ! Checks
                if ((tscr(1) /= 0.0_R8) .or. (tscr(size(tscr)) /= size(allc(i)%x)-1)) then 
                    ! Unexpected
                    call gdErrorHandler('SplitTMTubesTA: '  // &
                        'contour does not seem to start and end in a tube face, ' // &
                        'unexpected since tracing should start from face')
                end if 

                ! Construct the cell polygon
                tubecell = topomesh%tube%GetCell(tubeID(i))
                call ConstructTopologicalMeshCellPolygon(topomesh, &
                    tubecell(1), xcrda, ycrda)
                call cellpol%Construct(xcrda%Get(), ycrda%Get())

                ! Loop over all segments
                cc = 0
                wasfound = .false. 
                do j = 1, size(tscr)-1
                    ! Check if this segment should be considered
                    if (((tfaceind(j) == tubef(1)) .and. (tfaceind(j+1) == tubef(2))) .or. &
                        ((tfaceind(j) == tubef(2)) .and. (tfaceind(j+1) == tubef(1)))) then 
                            
                        ! Update the counter
                        cc = cc + 1

                        ! Extract the segment coordinates without end points
                        tempind = [(k, k = 1, size(allc(i)%x))]
                        tempind = pack(tempind, (tempind > tscr(j)+1.0_R8) .and. &
                            (tempind < tscr(j+1)+1.0_R8))
                        tempx = [allc(i)%x(tempind)]
                        tempy = [allc(i)%y(tempind)]

                        ! Check sizes
                        if (size(tempx) == 0) then 
                            ! Need to refine a bit
                            print *, 'SplitTMTubesTA: closed contour segment ' // & 
                                'does not have any vertices in between, attempting ' // & 
                                'to refine (results may be inaccurate...)'
                            tempx = real([(k, k = 0, 3)], kind=R8)/3.0_R8*(xint(j+1) - xint(j)) + xint(j)
                            tempy = real([(k, k = 0, 3)], kind=R8)/3.0_R8*(yint(j+1) - yint(j)) + yint(j)
                            tempx = tempx(2:size(tempx)-1)
                            tempy = tempy(2:size(tempy)-1)
                        end if 

                        ! Check if there are any vertices in the cell polygon
                        call Inpolygon(cellpol, tempx, tempy, isincellpol)
                        if (any(isincellpol)) then 
                            ! Check if already found a segment
                            if (wasfound) then 
                                call gdErrorHandler('SplitTMTubesTA: ' // & 
                                    'both segments appear to be located ' // &
                                    'inside of tube, unexpected')
                            end if 
                            wasfound = .true.
                            newx = [xint(j), tempx, xint(j+1)] 
                            newy = [yint(j), tempy, yint(j+1)]
                        end if 
                    end if
                end do 

                ! Sanity check: should've found two segments
                if (cc /= 2) then 
                    call gdErrorHandler('SplitTMTubesTA: closed contour ' // & 
                        'did not have two segments as expected, this may be a bug')
                end if 

                ! Sanity check: should've found a segment that lies within
                ! the polygon cell
                if (.not. wasfound) then 
                    call gdErrorHandler('SplitTMTubesTA: could not find ' // & 
                        'a closed contour segment that lies within the tube')
                end if 

                ! Keep only part inbetween intersections
                allc(i)%x = newx
                allc(i)%y = newy

                ! Keep this contour
                keepc(i) = .true.
                    
            else ! open contour
                ! Normally, there should be two contour parts for the 
                ! same tube, of which only one holds a segment that
                ! lies within the tube (i.e. of which there is a segment
                ! that starts at the starting face and ends at the tube 
                ! end face). Other contours are marked for deletion. 
                ! Check which intersection is the last intersection with
                ! the starting face (should be first one)

                ! Keep only intersections that intersect in one of the 
                ! tube's bounding radial faces
                allocate(keepf(size(tfaceind)))
                keepf = .false. 
                where (tfaceind == tubef(1) .or. tfaceind == tubef(2)) keepf = .true. 
                tscr = pack(tscr, keepf)
                xint = pack(xint, keepf)
                yint = pack(yint, keepf)
                tfaceind = pack(tfaceind, keepf)
                deallocate(keepf)

                ! Find start and end index that are subsequent. If none
                ! found, then this is not the contour we're looking for 
                ! and it should be deleted. 
                intersectind = findloc1D(tfaceind, tubef(1:2))

                ! Check if part is present
                if (intersectind == 0) then 
                    ! Remove contour and exit
                    keepc(i) = .false. 
                else
                    ! Keep only this segment of the contour
                    keepc(i) = .true.
                    allc(i)%startsaddle = 0
                    if (tscr(intersectind) == 0.0_R8) then 
                        call DeleteCurveSegment(allc(i)%x, allc(i)%y, &
                            [tscr(intersectind:intersectind+1)], 'both', [distfrac, -distfrac], .true., .false.)
                    else
                        call DeleteCurveSegment(allc(i)%x, allc(i)%y, &
                            [tscr(intersectind:intersectind+1)], 'both', [-distfrac, -distfrac], .false., .false.)
                    end if 
                end if 
                if (any(.not. isstartface(1:intersectind))) then 
                    ! Issue warning - we got intersections of the contour
                    ! with non-neighbour faces inbetween - this is unexpected
                    print *, 'InsertAlignedVesselParts: intersections ' // & 
                        'found for tangency point ', facevert(cface%Get(i)), &
                        'that occur inbetween intersections with contour and ' // & 
                        'starting boundary face. Unexpected, intersections ' // & 
                        'are removed'
                end if 
            
            end if 

            ! Housekeeping
            deallocate(isstartface)

        end do 

        ! Check if for all tubes contours were found. If not, issue 
        ! message (not per se an error since the topomesh will be fine)
        allocate(remtubeID(count(keepc)))
        remtubeID = pack(tubeID, keepc)
        call Setdiff(tubeID, remtubeID, difftubeID)
        if (size(difftubeID) /= 0) then 
            print *, 'InsertAlignedVesselParts: contours were not inserted ' // & 
                'for tubes: ', difftubeID 
        end if 
        deallocate(remtubeID)

        ! Housekeeping
        end associate

        ! Insert contours
        !================
        ! Insert
        call CleanContours(allc)
        newfsIDs = fsIDs%Get()
        if (options%writedebugoutput) then 
            allocate(allpc(ntpc))
            j = 0
            do i = 1, size(allc)
                !if (keepc(i)) then 
                    j = j + 1
                    call allpc(j)%Construct(allc(i)%x, allc(i)%y)
                !end if 
            end do 
            call tempps%Construct(allpc)
            call tempps%WriteData('avpcontours_all_beforeinsertion')
        end if 
        do i = 1, size(allc)
            if (keepc(i)) then 
                call InsertTopologicalMeshContour(topomesh, magneticField, &
                    allc(i), curvetypes%Get(i), newfsIDs(i))
            end if 
        end do
        topomesh%nFs = nfs

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp1', .false.)
        end if 

        ! Trim - normally, original tangency points etc shouldn't be removed
        ! so we can map the vertexmark logical index back
        call TrimTopologicalMesh(topomesh, magneticField, vessel)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp2', .false.)
        end if 

        ! Split boundaries (should probably not happen but ok)
        call SplitTopologicalMeshFaces(topomesh)

        ! Simplify
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Recompute marked vertex logical to deal with tubecase 2 later
        deallocate(vertexmark)
        allocate(vertexmark(topomesh%vert%ntot))
        vertexmark = .false.
        vertexmark(vertexmarkIDs) = .true.

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp3', .false.)
        end if 

        ! Recompute tubes
        !================
        
        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Tubes
        !======
        ! Flux tubes
        call AddTopologicalMeshTubes(topomesh, [TMfaceradID, TMfacebndID])

        ! Additional tube interconnection data
        call AddTopologicalMeshTubeData(topomesh)

        ! Adjust types of aligned boundaries 
        allocate(markedtpIDs(count(vertexmark)))
        markedtpIDs = pack([(k, k = 1, size(vertexmark))], vertexmark)
        do i = 1, topomesh%face%ntot
            ! Check if it contains any marked tangency points and 
            ! new flux surface IDs
            if (topomesh%face%type(i) == TMfacebndID) then 
                if (any(topomesh%face%vert(i, 1) == markedtpIDs) .and. &
                    any(topomesh%vert%fsID(topomesh%face%vert(i, 2)) == newfsIDs)) then 
                    ! Hedge for vertices that are remnants of out of vessel 
                    ! boundaries - have new flux surface ID but only two 
                    ! boundary faces as neighbours, so will be removed 
                    ! during simplification step
                    if (any( &
                        (topomesh%face%vert(i, 2) == topomesh%face%vert(:, 1) .or. &
                        topomesh%face%vert(i, 2) == topomesh%face%vert(:, 2)) .and. &
                        topomesh%face%type /= TMfacebndID)) then 
                        ! Change type
                        topomesh%face%type(i) = TMfacealbndID

                        ! Add ID
                        topomesh%face%fsID(i) = topomesh%vert%fsID(topomesh%face%vert(i, 1))
                    end if 
                elseif (any(topomesh%face%vert(i, 2) == markedtpIDs) .and. &
                        any(topomesh%vert%fsID(topomesh%face%vert(i, 1)) == newfsIDs)) then 
                        if (any( &
                            (topomesh%face%vert(i, 2) == topomesh%face%vert(:, 1) .or. &
                            topomesh%face%vert(i, 2) == topomesh%face%vert(:, 2)) .and. &
                        topomesh%face%type /= TMfacebndID)) then 
                        ! Change type
                        topomesh%face%type(i) = TMfacealbndID

                        ! Add ID
                        topomesh%face%fsID(i) = topomesh%vert%fsID(topomesh%face%vert(i, 1))
                    end if 
                end if 
            end if
        end do 

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp4', .false.)
        end if 

        ! Mark faces for deletion
        !========================
        allocate(remf(topomesh%face%ntot))
        remf = .false. 

        ! Loop a first time over all tubes
        do i = 1, topomesh%tube%ntot 
            ! Get tube faces
            tf = GetTMTubeFace(topomesh%tube, i)

            ! Check
            if (any(topomesh%face%type(tf) == TMfacealbndID)) then 
                ! Get neighbouring faces
                tfnb1 = GetTMTubeBndFace(topomesh%tube, i, 1)
                tfnb2 = GetTMTubeBndFace(topomesh%tube, i, 2)

                ! Check if first and/or second boundary has the 
                ! new flux surface ID - if so, delete the other
                hasnewfsIDnb1 = .false.
                hasnewfsIDnb2 = .false. 

                do j = 1, size(tfnb1)
                    if (any(topomesh%face%fsID(tfnb1(j)) == newfsIDs)) then 
                        hasnewfsIDnb1 = .true. 
                        exit ! no need to check further 
                    end if 
                end do
                do j = 1, size(tfnb2)
                    if (any(topomesh%face%fsID(tfnb2(j)) == newfsIDs)) then 
                        hasnewfsIDnb2 = .true. 
                        exit ! no need to check further 
                    end if 
                end do

                ! Sanity checks
                if (.not. hasnewfsIDnb1 .and. .not. hasnewfsIDnb2) then 
                    ! This can only happen if the tube is fully aligned, 
                    ! which shouldn't be the case here anymore

                    tfnbv1 = [topomesh%face%vert(tfnb1, 1), topomesh%face%vert(tfnb1, 2)]
                    tfnbv2 = [topomesh%face%vert(tfnb2, 1), topomesh%face%vert(tfnb2, 2)]
                    print *, 'tube: ', i
                    print *, 'vertices: ', tfnbv1, tfnbv2
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // & 
                        'tube has aligned boundary face but no new flux ' // & 
                        'surface as boundary - unexpected.')
                elseif (hasnewfsIDnb1 .and. hasnewfsIDnb2) then 
                    print *, 'tube: ', i
                    call WriteTopologicalMesh(topomesh, 'topomesh_error')
                    call gdErrorHandler('InsertAlignedVesselParts: ' // & 
                        'tube has aligned boundary face but both boundaries have new flux ' // & 
                        'surface as boundary - unexpected.')
                elseif (hasnewfsIDnb1) then 
                    ! Delete the second faces
                    remf(tfnb2) = .true.
                else 
                    ! Delete the first faces
                    remf(tfnb1) = .true. 
                end if 
            end if 
        end do 

        ! Adjust some boundaries again to be aligned based on which 
        ! faces that will be removed
        do i = 1, topomesh%vert%ntot
            ! Check if the vertex has the same flux surface ID as one 
            ! of the marked tangency points and if it's a boundary vertex
            if (any(topomesh%vert%fsID(i) == topomesh%vert%fsID(markedtpIDs)) .and. &
                (topomesh%vert%type(i) == TMvertexbndID)) then 
                ! Get the vertex face neighbours
                tvf = topomesh%vert%GetFace(i)

                ! If it has an aligned boundary, and if one of its
                ! neighbouring faces is an aligned face that 
                ! was marked for deletion, then  all other boundary faces 
                ! should become aligned too - otherwise non-conforming
                ! tubes may exist
                if (any(topomesh%face%type(tvf) == TMfacealbndID) .and. &
                    any(remf(tvf))) then 
                    where (topomesh%face%type(tvf) == TMfacebndID) 
                        topomesh%face%type(tvf) = TMfacealbndID
                        topomesh%face%fsID(tvf) = topomesh%vert%fsID(i)
                    end where
                end if 
            end if 
        end do 

        ! Check if we need to remove vertices (only if they have no faces left)
        allocate(remv(topomesh%vert%ntot))
        remv = .false. 
        do i = 1, topomesh%vert%ntot
            tvf = topomesh%vert%GetFace(i)
            if (all(remf(tvf))) then 
                remv(i) = .true. 
            end if 
        end do 

        ! Remove faces
        call RemoveTopologicalMeshFaceLogical(topomesh, remf)

        ! Remove vertices
        call RemoveTopologicalMeshVertexLogical(topomesh, remv)

        ! Simplify 
        call SimplifyTopologicalMeshFaces(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_duringavp5', .false.)
        end if 

        ! Recompute interconnection data
        !===============================
        ! Vertex faces
        call AddTopologicalMeshVertexFaces(topomesh)

        ! Data 
        call AddTopologicalMeshData(topomesh)

        ! Add cells
        call AddTopologicalMeshCells(topomesh)

        ! Recompute interconnection data
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Do temporary writing
        if (options%writedebugoutput) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_afteravp')
        end if 
    
    end subroutine

    ! Wrapper for tube pair getter
    subroutine GetMergeTubePairsWrapper(tmadaptor, topomesh, hftubes, &
        lftubes, includetube)

        ! Description
        !============
        ! A simple wrapper to determine in which way tube pairs should
        ! be constructed (not too much choice, but still). Currently, 
        ! only two options supported: 'extensive' and 'local'. The 
        ! extensive option explores all neighbouring tubes and merges
        ! over a much larger surface than the local mode. This typically
        ! leads to more tubes remaining in the topomesh. The local mode
        ! only looks at two tubes at the same time, and can therefore 
        ! merge more tubes than the extensive method. However, this 
        ! may lead to possibly undesired side effects and 'uglier' 
        ! tubes. See the dedicated subroutines for additional explanation

        ! Arguments
        class(TopomeshAdaptorUDT)           :: tmadaptor
        class(TopomeshUDT), intent(in)      :: topomesh
        type(IntegerDynamicArrayUDT), allocatable, intent(out)  :: &
            hftubes(:), lftubes(:)
        logical, dimension(:), intent(in), optional     :: includetube

        ! Initialize
        !===========
        select case (tmadaptor%mergetubemeth)

        case ('extensive')

            if (present(includetube)) then 
                call tmadaptor%GetMergeTubePairsExtensive(topomesh, hftubes, &
                    lftubes, includetube)
            else
                call tmadaptor%GetMergeTubePairsExtensive(topomesh, hftubes, &
                    lftubes)
            end if
            
        case ('local')

            if (present(includetube)) then 
                call tmadaptor%GetMergeTubePairsLocal(topomesh, hftubes, &
                    lftubes, includetube)
            else
                call tmadaptor%GetMergeTubePairsLocal(topomesh, hftubes, &
                    lftubes)
            end if
            
        case default 

            call gdErrorHandler('GetMergeTubePairsWrapper: unknown tube ' // & 
                'merging method')

        end select

    end subroutine

    ! Extensive tube pair getter for merging
    subroutine GetMergeTubePairsExtensiveTA(tmadaptor, topomesh, hftubes, lftubes, &
        includetube)

        ! Description
        !============
        ! This routine returns sets of tubes (resp. high and low flux
        ! tubes) that are connected to each other on the low flux 
        ! boundary of the high flux neighbours and vice versa. Note that
        ! the search is exhaustive, in the sense that for each found 
        ! tube on e.g. the low side, all high flux neighbours are taken
        ! and again their low flux neighbours are taken. This process is
        ! repeated until no new tubes are found. Note that by definition 
        ! the hftubes cannot be empty if the lftubes are nonempty and 
        ! vice versa. Since this algorithm essentially assumes an acyclic
        ! graph, we need to hedge separately for cycles occuring in the 
        ! topomesh due to vessel structures. To this end, we identify 
        ! the cycles using a dedicated routine and afterwards check which
        ! tube pairs may be considered for merging. This should prevent
        ! the occurrence of detrimental cycle topologies that cannot 
        ! adequately be merged anymore (as long as they were not 
        ! initially present, which shouldn't be the case if the 
        ! topological mesh generation step went well)

        ! Note: optionally, the logical 'includetube' can be passed, 
        ! which should be a 1-by-ntubes logical indicating which tubes
        ! to consider. If not passed, all possible pairs are computed. 

        ! Algorithm
        !----------
        ! 1)    Take a tube
        ! 2)    Determine all high flux neighbours of the tube. For each 
        !       high flux neighbour, take its low flux neighbours and add.
        !       Repeat this step for each low flux neighbour, until all
        !       tubes have been reached. This loop can only be finite if
        !       there are no simple cycles present (this is checked for 
        !       and an error will be throw. However, it should be prevented
        !       by the explicit cycle checking done afterwards)
        ! 3)    Check if a set of tubes on high and low flux side could 
        !       be found. If this is the case, then proceed to 4). Otherwise,
        !       no pair could be found -> go to 5)
        ! 4)    Hedge for cycles: for each tube in the pair that is part
        !       of a cycle (or multiple cycles), check if any neighbour
        !       is a bounding tube. If this is the case, the pair can 
        !       only be considered if one of the following conditions hold:
        !       for each tube of the pair that is in the cycle, all 
        !       neighbouring tubes that belong to that cycle should 
        !       either be non-bounding tubes or bounding tubes that are
        !       already included in the pair
        !       OR
        !       each non-bounding tube of the cycle should have only
        !       bounding tubes as high and low field neighbours OR 
        !       neighbours that are not part of the cycle (otherwise, 
        !       merging may connect two bounding tubes directly, which 
        !       will prevent other merging operations to take place and 
        !       should hence be avoided).
        !       If this is not the case, we should be able to proceed
        ! 5)    Repeat 3)-4), but start with the low flux neighbours and 
        !       determine the high flux neighbours.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)           :: tmadaptor
        class(TopomeshUDT), intent(in)      :: topomesh
        type(IntegerDynamicArrayUDT), allocatable, intent(out)  :: &
            hftubes(:), lftubes(:)
        logical, dimension(:), intent(in), optional     :: includetube

        ! Auxiliary
        integer(I8)                             :: tpc, ind
        integer(I8), allocatable, dimension(:)  :: tnb, tv, thft, tlft, &
            ct, ctind, tct, tbtind, thfnb, tlfnb, thfnba, tlfnba, tnba
        logical                                 :: addpair
        logical, allocatable, dimension(:)      :: ishfnb, islfnb, &
            tracetubehf, tracetubelf, istubefound, ismarked, isincycle, &
            isboundingtube, isactivetube, iscycleboundingtube, &
            keeptubepair
        type(IntegerDynamicArrayUDT)            :: hft, lft
        type(IntegerDynamicArrayUDT), allocatable   :: tube1ida(:), &
            tube2ida(:), cycletubes(:), cycleboundingtubeind(:)
        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Unpack for ease at this stage
        associate(&
            face    => topomesh%face,   &
            tube    => topomesh%tube    &
            )

        ! Check optional arguments
        if (present(includetube)) then 
            ismarked = includetube
        else
            allocate(ismarked(tube%ntot))
            ismarked = .true.
        end if 

        ! Determine cycles
        !=================
        ! Initialize
        allocate(isincycle(tube%ntot), isboundingtube(tube%ntot))
        isincycle = .false.
        isboundingtube = .false. 

        ! Call dedicated routine
        call tmadaptor%GetMergeTubeCycles(topomesh, cycletubes, cycleboundingtubeind)

        ! Process information for easier use later on
        do i = 1, size(cycletubes)
            ! Get the current tubes and indices
            ct = cycletubes(i)%Get()
            ctind = cycleboundingtubeind(i)%Get()

            ! Set logicals
            isincycle(ct) = .true.
            isboundingtube(ct(ctind)) = .true. 
        end do

        ! Determine mergeable tube pairs
        !===============================
        ! Construct tube pairs
        tpc = 0 ! tube pair counter
        allocate(ishfnb(tube%ntot), islfnb(tube%ntot), tracetubehf(tube%ntot), &
            tracetubelf(tube%ntot), istubefound(tube%ntot))
        allocate(tube1ida(2*count(ismarked)), tube2ida(2*count(ismarked))) ! possibly too big
        ishfnb = .false. 
        islfnb = .false. 
        do i = 1, tube%ntot
            ! Skip in case tube was not makred
            if (.not. ismarked(i)) then 
                cycle
            end if

            ! High flux side
            !---------------
            ! Check if it's not already a low flux tube neighbour in 
            ! another set - in that case it is already part of a tubeset
            ! as low flux tube
            if (.not. islfnb(i)) then 
                ! Initialize
                tracetubelf = .false. 
                tracetubehf = .false. 
                istubefound = .false. 
                hft = ConstructIntegerDynamicArray()
                lft = ConstructIntegerDynamicArray()
            
                ! Set as found for low flux neighbour
                islfnb(i) = .true. 
                istubefound(i) = .true.

                ! Mark as tube to trace high flux neighbours
                tracetubehf(i) = .true. 

                ! Add to lft
                call lft%Append(i)

                ! Find all other tubes 
                do while (any(tracetubehf) .or. any(tracetubelf))
                    ! Trace high flux side
                    ind = findloc(tracetubehf, .true., 1, back=.false.)
                    if (ind /= 0) then 
                        ! Remove from tracing
                        tracetubehf(ind) = .false. 

                        ! Find tubes at high flux side of this tube
                        tnb = tube%GetHighFluxNeig(ind)

                        ! Retain only tubes that were not found already
                        tnb = pack(tnb, .not. istubefound(tnb))

                        ! Check if none of these tubes were high flux 
                        ! tubes already - normally this should not be the case
                        if (any(ishfnb(tnb))) then 
                            tv = [tube%GetBndVert(i, 1), tube%GetBndVert(i, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            call WriteTopologicalMesh(topomesh, 'topomesh_error')
                            call gdErrorHandler('GetMergeTubePairs: ' // & 
                                'tubes found that are already marked as ' // & 
                                'high flux tubes, unexpected')
                        end if 

                        ! Mark tubes to trace their low flux side
                        tracetubelf(tnb) = .true. 

                        ! Mark tubes as found
                        istubefound(tnb) = .true.
                        ishfnb(tnb) = .true.

                        ! Add to hft
                        call hft%Append(tnb)
                    end if 

                    ! Trace low flux side
                    ind = findloc(tracetubelf, .true., 1, back=.false.)
                    if (ind /= 0) then 
                        ! Remove from tracing
                        tracetubelf(ind) = .false. 

                        ! Find tubes at high flux side of this tube
                        tnb = tube%GetLowFluxNeig(ind)

                        ! Retain only tubes that were not found already
                        tnb = pack(tnb, .not. istubefound(tnb))

                        ! Check if none of these tubes were high flux 
                        ! tubes already - normally this should not be the case
                        if (any(islfnb(tnb))) then 
                            tv = [tube%GetBndVert(i, 1), tube%GetBndVert(i, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            call WriteTopologicalMesh(topomesh, 'topomesh_error.dat')
                            call gdErrorHandler('GetMergeTubePairs: ' // & 
                                'tubes found that are already marked as ' // & 
                                'low flux tubes, unexpected')
                        end if 

                        ! Mark tubes to trace their high flux side
                        tracetubehf(tnb) = .true. 
                        
                        ! Mark tubes as found
                        istubefound(tnb) = .true.
                        islfnb(tnb) = .true.

                        ! Add to lft
                        call lft%Append(tnb)
                    end if 
                end do 

                ! Add pair if high flux neighbours were found and if 
                ! none of the high flux neighbours has another high 
                ! flux neighbour as high flux neighbour 
                if (hft%Size() /= 0) then 
                    ! Check high flux neighbours of high flux neighbours
                    thft = hft%Get()
                    addpair = .true. 
                    !do j = 1, size(thft)
                    !    tnb = tube%GetHighFluxNeig(thft(j))
                    !    if (size(GetCommonElements(tnb, thft)) /= 0) then
                    !        addpair = .false.
                    !        exit
                    !    end if 
                    !end do  

                    ! Add pair if possible
                    if (addpair) then 
                        ! Increment tube pair counter
                        tpc = tpc + 1

                        ! Initialize 
                        tube1ida(tpc) = ConstructIntegerDynamicArray(hft%Get())
                        tube2ida(tpc) = ConstructIntegerDynamicArray(lft%Get())
                    else
                        ! Allow pairs to be found again
                        islfnb(lft%Get()) = .false.
                        ishfnb(hft%Get()) = .false.
                        islfnb(i) = .false.
                    end if 
                else
                    ! Allow lft tubes to be found again
                    islfnb(lft%Get()) = .false.
                    islfnb(i) = .false.
                end if 
            end if 

            ! Low flux side
            !--------------
            ! Check if it's not already a high flux tube neighbour in 
            ! another set - in that case it is already part of a tubeset
            ! as low flux tube
            if (.not. ishfnb(i)) then 
                ! Initialize
                tracetubelf = .false. 
                tracetubehf = .false. 
                istubefound = .false. 
                hft = ConstructIntegerDynamicArray()
                lft = ConstructIntegerDynamicArray()
            
                ! Set as found for high flux neighbour
                ishfnb(i) = .true. 
                istubefound(i) = .true.

                ! Mark as tube to trace high flux neighbours
                tracetubelf(i) = .true. 

                ! Add to hft
                call hft%Append(i)

                ! Find all other tubes 
                do while (any(tracetubehf) .or. any(tracetubelf))
                    ! Trace high flux side
                    ind = findloc(tracetubehf, .true., 1, back=.false.)
                    if (ind /= 0) then 
                        ! Remove from tracing
                        tracetubehf(ind) = .false. 

                        ! Find tubes at high flux side of this tube
                        tnb = tube%GetHighFluxNeig(ind)

                        ! Retain only tubes that were not found already
                        tnb = pack(tnb, .not. istubefound(tnb))

                        ! Check if none of these tubes were high flux 
                        ! tubes already - normally this should not be the case
                        if (any(ishfnb(tnb))) then 
                            tv = [tube%GetBndVert(i, 1), tube%GetBndVert(i, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            call WriteTopologicalMesh(topomesh, 'topomesh_error.dat')
                            call gdErrorHandler('GetMergeTubePairs: ' // & 
                                'tubes found that are already marked as ' // & 
                                'high flux tubes, unexpected')
                        end if 

                        ! Mark tubes to trace their low flux side
                        tracetubelf(tnb) = .true. 

                        ! Mark tubes as found
                        istubefound(tnb) = .true.
                        ishfnb(tnb) = .true.

                        ! Add to hft
                        call hft%Append(tnb)
                    end if 

                    ! Trace low flux side
                    ind = findloc(tracetubelf, .true., 1, back=.false.)
                    if (ind /= 0) then 
                        ! Remove from tracing
                        tracetubelf(ind) = .false. 

                        ! Find tubes at high flux side of this tube
                        tnb = tube%GetLowFluxNeig(ind)

                        ! Retain only tubes that were not found already
                        tnb = pack(tnb, .not. istubefound(tnb))

                        ! Check if none of these tubes were high flux 
                        ! tubes already - normally this should not be the case
                        if (any(islfnb(tnb))) then 
                            tv = [tube%GetBndVert(i, 1), tube%GetBndVert(i, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            tv = [tube%GetBndVert(5, 1), tube%GetBndVert(5, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            tv = [tube%GetBndVert(11, 1), tube%GetBndVert(11, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            tv = [tube%GetBndVert(20, 1), tube%GetBndVert(20, 2)]
                            print *, 'tube vertices of current tube: ', tv
                            call WriteTopologicalMesh(topomesh, 'topomesh_error.dat')
                            call gdErrorHandler('GetMergeTubePairs: ' // & 
                                'tubes found that are already marked as ' // & 
                                'low flux tubes, unexpected')
                        end if 

                        ! Mark tubes to trace their high flux side
                        tracetubehf(tnb) = .true. 
                        
                        ! Mark tubes as found
                        istubefound(tnb) = .true.
                        islfnb(tnb) = .true.

                        ! Add to lft
                        call lft%Append(tnb)
                    end if 
                end do 

                ! Add pair if low flux neighbours were found
                if (lft%Size() /= 0) then 
                    ! Check low flux neighbours of low flux neighbours
                    addpair = .true.
                    tlft = lft%Get()
                    !do j = 1, size(tlft)
                    !    tnb = tube%GetLowFluxNeig(tlft(j))
                    !    if (size(GetCommonElements(tnb, tlft)) /= 0) then 
                    !        addpair = .false.
                    !        exit
                    !    end if 
                    !end do 
                    
                    ! Check if we can add the pair
                    if (addpair) then 
                        ! Increment tube pair counter
                        tpc = tpc + 1

                        ! Initialize 
                        tube1ida(tpc) = ConstructIntegerDynamicArray(hft%Get())
                        tube2ida(tpc) = ConstructIntegerDynamicArray(lft%Get())
                    else

                        ! Allow pairs to be found again
                        islfnb(lft%Get()) = .false.
                        ishfnb(hft%Get()) = .false.
                        ishfnb(i) = .false.
                    end if 
                else
                    ! Allow lft tubes to be found again
                    islfnb(hft%Get()) = .false.
                    ishfnb(i) = .false.
                end if 
            end if 
        end do 

        ! If no pairs were found, exit
        if (tpc == 0) then 
            allocate(hftubes(0), lftubes(0))
            return 
        end if 

        ! Remove pairs based on cycles
        allocate(isactivetube(tube%ntot), keeptubepair(tpc))
        keeptubepair = .true. 
        do i = 1, tpc
            ! Unpack
            thft = tube1ida(i)%Get()
            tlft = tube2ida(i)%Get()

            ! Check if there are any tubes of a cycle in the pair
            if (.not. any(isincycle([tlft, thft]))) then 
                ! Keep this pair, skip
                cycle
            end if 

            ! Check if there are any bounding tubes present
            if (.not. any(isboundingtube([tlft, thft]))) then 
                ! Keep this pair, skip
                cycle
            end if 

            ! Bounding tubes are present, so we need to do thorough 
            ! checks for each cycle
            isactivetube = .false.
            isactivetube([thft, tlft]) = .true. 
            do j = 1, size(cycletubes)
                ! Check if the tube pair was already deleted, if so, exit the loop
                if (.not. keeptubepair(i)) then 
                    exit
                end if

                ! Get the current active cycle tubes
                tct = cycletubes(j)%Get()
                tbtind = cycleboundingtubeind(j)%Get()

                ! Check if any cycle tubes are in the current tube pair
                if (.not. any(isactivetube(tct))) then 
                    ! Skip
                    cycle
                end if 

                ! Check if any active cycle tubes are bounding tubes
                if (.not. any(isactivetube(tct(tbtind)))) then 
                    ! Skip
                    cycle
                end if

                ! If we got here, then we need to check the conditions 
                ! in step 4) of the algorithm
                allocate(iscycleboundingtube(tube%ntot))
                iscycleboundingtube = .false.
                iscycleboundingtube(tct(tbtind)) = .true.
                
                ! Get all high/low flux neighbours of non-bounding tubes 
                ! of this cycle
                allocate(thfnb(0), tlfnb(0), thfnba(0), tlfnba(0))
                do k = 1, size(tct)
                    ! Skip if it is a bounding tube - these don't have 
                    ! to be checked
                    if (iscycleboundingtube(tct(k))) then 
                        cycle
                    end if 
                    thfnb = [thfnb, tube%GetHighFluxNeig(tct(k))]
                    tlfnb = [tlfnb, tube%GetLowFluxNeig(tct(k))]
                    if (isactivetube(tct(k))) then 
                        thfnba = [thfnba, tube%GetHighFluxNeig(tct(k))]
                        tlfnba = [tlfnba, tube%GetLowFluxNeig(tct(k))]
                    end if 


                    ! Sanity check
                    if (size([thfnb, tlfnb]) < 1) then 
                        ! Apparently cycle tube has no neighbours, this 
                        ! should not be possible
                        print *, 'tube vertices: ', tube%GetBndVert(tct(k), 1), tube%GetBndVert(tct(k), 2)
                        call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                        call gdErrorHandler('GetMergeTubePairs: cycle tube ' // &
                            'does not have any cycle neighbours, unexpected')
                    end if 
                end do 

                ! Keep only neighbours that are part of the current 
                ! cycle
                thfnb = GetCommonElements(thfnb, tct)
                tlfnb = GetCommonElements(tlfnb, tct)
                tnb = [thfnb, tlfnb]
                thfnba = GetCommonElements(thfnba, tct)
                tlfnba = GetCommonElements(tlfnba, tct)
                tnba = [thfnba, tlfnba]

                ! Check second condition: all neighbours are tube neighbours
                if (all(iscycleboundingtube(tnb)) .and. all(isactivetube(tct))) then 
                    ! Keep pair
                else
                    ! first condition: all neighbours of active cycle 
                    ! tubes should be considered
                    do k = 1, size(tnba)
                        if (iscycleboundingtube(tnba(k)) .and. .not. isactivetube(tnba(k))) then 
                            keeptubepair(i) = .false.
                            exit
                        end if 
                    end do 
                end if 

                ! Housekeeping
                deallocate(iscycleboundingtube, thfnb, tlfnb, thfnba, tlfnba)
            end do 
        end do

        ! Set output
        hftubes = tube1ida(1:tpc)
        lftubes = tube2ida(1:tpc)
        if (.not. all(keeptubepair)) then 
            hftubes = pack(hftubes, keeptubepair)
            lftubes = pack(lftubes, keeptubepair)
        end if 

        ! Housekeeping
        end associate

    end subroutine

    ! Local tube pair getter for merging
    subroutine GetMergeTubePairsLocalTA(tmadaptor, topomesh, hftubes, lftubes, &
        includetube)

        ! Description
        !============
        ! This routine returns sets of tubes (resp. high and low flux
        ! tubes) that are connected to each other on the low flux 
        ! boundary of the high flux neighbours and vice versa. In this
        ! routine, we check for each tube to be considered the high and 
        ! low flux neighbours. Each tube pair formed that way then 
        ! simply becomes a tube pair (so in principle, all hftubes and 
        ! lftubes will always have one tube, but to be consistent with 
        ! previous implementation we still return all as an array of 
        ! integer dynamic arrays). We do exclude pairs that would consist
        ! of only bounding tubes in a cycle to prevent merging away 
        ! 'internal' vessel structures. This approach will likely lead
        ! to more aligned vessel faces in the topomesh (which may 
        ! actually be desired from the gridding point of view)
        ! Note: optionally, the logical 'includetube' can be passed, 
        ! which should be a 1-by-ntubes logical indicating which tubes
        ! to consider. If not passed, all possible pairs are computed. 

        ! Algorithm
        !----------
        ! 1)    Construct initial tube pairs based on the graph of the
        !       topomesh - each graph edge is considered to be a pair
        ! 2)    For each pair, hedge for cycles by checking if both 
        !       pairs are bounding tubes in the same cycle. This kind of
        !       merge should only be possible if the cycle only contains
        !       bounding tubes (otherwise, the pair should be deleted)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)           :: tmadaptor
        class(TopomeshUDT), intent(in)      :: topomesh
        type(IntegerDynamicArrayUDT), allocatable, intent(out)  :: &
            hftubes(:), lftubes(:)
        logical, dimension(:), intent(in), optional     :: includetube

        ! Auxiliary
        integer(I8)                             :: thft, tlft, tpc
        integer(I8), allocatable, dimension(:)  :: ct, ctind, tct, &
            tbtind, tube1, tube2
        logical, allocatable, dimension(:)      :: ismarked, isincycle, &
            isboundingtube, isactivetube, keeptubepair, iscurrentboundingtube
        type(IntegerDynamicArrayUDT), allocatable   :: cycletubes(:), &
            cycleboundingtubeind(:)
        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Unpack for ease at this stage
        associate(&
            face    => topomesh%face,   &
            tube    => topomesh%tube    &
            )

        ! Check optional arguments
        if (present(includetube)) then 
            ismarked = includetube
        else
            allocate(ismarked(tube%ntot))
            ismarked = .true.
        end if 

        ! Determine cycles
        !=================
        ! Initialize
        allocate(isincycle(tube%ntot), isboundingtube(tube%ntot))
        isincycle = .false.
        isboundingtube = .false. 

        ! Call dedicated routine
        call tmadaptor%GetMergeTubeCycles(topomesh, cycletubes, cycleboundingtubeind)

        ! Process information for easier use later on
        do i = 1, size(cycletubes)
            ! Get the current tubes and indices
            ct = cycletubes(i)%Get()
            ctind = cycleboundingtubeind(i)%Get()

            ! Set logicals
            isincycle(ct) = .true.
            isboundingtube(ct(ctind)) = .true. 
        end do

        ! Determine mergeable tube pairs
        !===============================
        tube1 = topomesh%tube%graph%v(topomesh%tube%graph%ev1)
        tube2 = topomesh%tube%graph%v(topomesh%tube%graph%ev2)

        ! Construct tube pairs
        tpc = size(tube1) ! tube pair counter

        ! If no pairs were found, exit
        if (tpc == 0) then 
            allocate(hftubes(0), lftubes(0))
            return 
        end if 

        ! Remove pairs based on cycles
        allocate(isactivetube(tube%ntot), keeptubepair(tpc), &
            iscurrentboundingtube(tube%ntot))
        keeptubepair = .true. 
        do i = 1, tpc
            ! Unpack
            thft = tube1(i)
            tlft = tube2(i)

            ! Check if there are any tubes of a cycle in the pair
            if (.not. any(isincycle([tlft, thft]))) then 
                ! Keep this pair, skip
                cycle
            end if 

            ! Check if there are any bounding tubes present
            if (.not. any(isboundingtube([tlft, thft]))) then 
                ! Keep this pair, skip
                cycle
            end if 

            ! Bounding tubes are present, so we need to do thorough 
            ! checks for each cycle
            isactivetube = .false.
            isactivetube([thft, tlft]) = .true. 
            do j = 1, size(cycletubes)
                ! Check if the tube pair was already deleted, if so, exit the loop
                if (.not. keeptubepair(i)) then 
                    exit
                end if

                ! Get the current active cycle tubes
                tct = cycletubes(j)%Get()
                tbtind = cycleboundingtubeind(j)%Get()

                ! Check if any cycle tubes are in the current tube pair
                if (.not. any(isactivetube(tct))) then 
                    ! Skip
                    cycle
                end if 

                ! Check if any active cycle tubes are bounding tubes
                if (.not. any(isactivetube(tct(tbtind)))) then 
                    ! Skip
                    cycle
                end if

                ! Check if all active tubes are bounding tubes
                iscurrentboundingtube = .false.
                iscurrentboundingtube(tct(tbtind)) = .true.
                if (.not. all(iscurrentboundingtube([thft, tlft]))) then 
                    cycle
                end if 

                ! Bounding tubes are present, so can only include pair 
                ! if all are bounding tubes in this cycle
                if (.not. (size(tbtind )== size(tct))) then 
                    ! Mark 
                    keeptubepair(i) = .false. 
                end if
            end do 
        end do

        ! Set output
        allocate(hftubes(count(keeptubepair)), lftubes(count(keeptubepair)))
        j = 0
        do i = 1, size(keeptubepair)
            if (keeptubepair(i)) then
                j = j + 1
                hftubes(j) = ConstructIntegerDynamicArray([tube1(i)])
                lftubes(j) = ConstructIntegerDynamicArray([tube2(i)])
            end if 
        end do 

        ! Housekeeping
        end associate

    end subroutine

    ! Tube cycle getter for merging
    subroutine GetTopomeshTubeCyclesTA(tmadaptor, topomesh, tubes, boundingtubeind)

        ! Description
        !============
        ! This is a *very* specific routine to determine cycles and 
        ! other relevant information for topological mesh adaptations. 
        ! We fully exploit the assumption of only having (nested) closed
        ! boundary structures and the basic properties of the magnetic 
        ! field (i.e. that there are low/high flux sides and that 
        ! hence regions that are not bounded by any structure cannot 
        ! lead to cycles). This should be faster/easier than an 
        ! algorithm for cycle detection in graphs, which would be more
        ! general. We return all tubes adjacent in this way to the 
        ! polygon. Additionally, we also return the tubes that are only
        ! adjacent by aligned vessel faces or type 2 tangency points, 
        ! since these are 'bounding' tubes of the polygon (these should 
        ! not occur together in a merge, unless only these two tubes are present, 
        ! in which case this isn't a cycle anymore). Note that the amount
        ! of bounding tubes can be larger than two, since X-points may
        ! be present inside the closed polygon. 

        ! Notes
        !======
        ! Note 1: the tubes are *not* returned in any sorted way

        ! Algorithm
        !==========
        ! 1)    Find all closed polygons formed by (aligned) boundary 
        !       polygons. Only tubes around these polygons can lead to 
        !       cycles in the topomesh
        ! 2)    Determine the closed polygon nestedness level based
        !       on data from the vessel. Only odd levels should be 
        !       taken, as these have the tubes at the 'outside' and 
        !       only these can lead to cycles in the graph. 
        ! 3)    Take a closed polygon and continue to 4). If none are 
        !       left, exit
        ! 4)    Find tubes that are adjacent to the faces and vertices
        !       of each closed polygon. 
        ! 5)    If only two tubes are present, or if no bounding tubes 
        !       exist, it is not a cycle. Go to 3). Otherwise, continue to 6)
        ! 6)    Mark these tubes as bounding tubes and add data to output.
        !       Go to 3) 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)           :: tmadaptor
        class(TopomeshUDT), intent(in)      :: topomesh
        type(IntegerDynamicArrayUDT), intent(out), allocatable  :: &
            boundingtubeind(:), tubes(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)      :: tf, sortind, tpf, &
            tpfv, tubev, tubef, polygonlevel
        integer(I8), allocatable, dimension(:, :)   :: tfv, temp
        logical, allocatable, dimension(:)          :: isbranchingpolygon, &
            ispolygonstart
        type(IntegerDynamicArrayUDT)                :: thistubes, thisboundaryind

        ! Loop
        integer(I8)                                 :: i, si, ei, nc

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            tube    => topomesh%tube,   &
            cell    => topomesh%cell,   &
            face    => topomesh%face,   &
            vert    => topomesh%vert    &
            )

        ! Find polygons
        !==============
        ! Get all topomesh vessel faces & vert
        tf = topomesh%GetVesselFaceIDs()
        tfv = topomesh%face%vert(tf, :)

        ! Sort faces
        allocate(sortind(size(tf)), ispolygonstart(size(tf)), &
            isbranchingpolygon(size(tf)))
        call SortPolygonEdges(tfv, size(tf), sortind, ispolygonstart, &
            isbranchingpolygon)
        tfv = tfv(sortind, :)
        tf = tf(sortind)

        ! Sanity checks
        if (any(isbranchingpolygon)) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_error')
            call gdErrorHandler('GetTopomeshTubeCycles: boundary faces ' // & 
                'of topomesh seem to form branching polygons, not supported')
        end if 

        ! Loop over all polygons
        !=======================
        ! Initialize
        si = 0
        ei = 0
        nc = 0
        allocate(tubes(count(ispolygonstart)), boundingtubeind(count(ispolygonstart))) ! overestimation
        do while (.true.)
            ! Take polygon
            !-------------
            ! Find the next polygon faces 
            si = findloc(ispolygonstart(ei+1:), .true., 1, back=.false.) + ei 
            if (si == ei) then 
                ! Reached the end - exit
                exit
            end if 
            if (si == size(ispolygonstart)) then 
                ei = si 
            else
                ei = findloc(ispolygonstart(si+1:), .true., 1, back=.false.) + si 
                if (ei == si) then 
                    ! Last polygon
                    ei = size(ispolygonstart)
                else
                    ei = ei - 1
                end if 
            end if 
            tpf = tf(si:ei)

            ! Check if the polygon is closed, otherwise skip
            call ExtractPolygonVertices(face%vert(tpf, :), size(tpf), &
                tpfv)
            if (tpfv(1) /= tpfv(size(tpfv))) then
                ! Not a closed polygon  
                cycle
            end if

            ! Check if the polygon level is odd (just evaluate at vertex
            ! locations)
            call tmadaptor%vessel%exactplfvessel%EvaluateLabel(topomesh%vert%x(tpfv), &
                topomesh%vert%y(tpfv), temp)
            call Unique(temp(:, 4), polygonlevel)
            if (size(polygonlevel) > 1) then 
                call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                call tmadaptor%vessel%exactplfvessel%VisualizeLabel('levelset_label_error_1', labelindin=1)
                call tmadaptor%vessel%exactplfvessel%VisualizeLabel('levelset_label_error_1', labelindin=4)
                call gdErrorHandler('GetTopomeshTubeCyclesTA: multiple polygon levels found, unexpected. ' // & 
                    'Check if polygon levelset function for vessel was correctly constructed, or ' // & 
                    'if vessel polygon parts are almost intersecting, which may cause this behavior.')
            end if
            if (mod(polygonlevel(1), 2) == 0) then 
                ! Even polygon, skip
                cycle
            end if  


            ! If we got here, we can initialize the temporary arrays
            thistubes = ConstructIntegerDynamicArray()
            thisboundaryind = ConstructIntegerDynamicArray()

            ! Find adjacent tubes
            !--------------------
            ! Simply check vertices
            do i = 1, tube%ntot
                ! Get tube vertices
                tubev = [tube%GetBndVert(i, 1), tube%GetBndVert(i, 2)]

                ! Check for common elements
                if (size(GetCommonElements(tubev, tpfv)) > 0) then 
                    ! Append
                    call thistubes%Append(i)

                    ! Check if the tube has any radial faces in common. 
                    ! If this is not the case, the tube is a bounding tube
                    tubef = tube%GetFace(i)
                    if (size(GetCommonElements(tubef, tpf)) == 0) then 
                        ! Append
                        call thisboundaryind%Append(thistubes%Size()) ! should be last tube
                    end if 
                end if 
            end do 

            ! Add
            !----
            ! Only if it is a cycle
            if (thistubes%Size() <= 2) then 
                cycle
            end if 
            if (thisboundaryind%Size() == 0) then 
                cycle
            end if 

            ! Add
            nc = nc + 1
            tubes(nc) = thistubes
            boundingtubeind(nc) = thisboundaryind

        end do 

        ! Trim output
        tubes = tubes(1:nc)
        boundingtubeind = boundingtubeind(1:nc)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                   EQUILIBRIUM CHARACTERIZATION                   !
    !------------------------------------------------------------------!

    ! 2D extrema 
    subroutine TraceExtrema2D(xe, ye, fe, typee, fieldtracer, &
        magneticField, donewton)

        ! Description
        !============
        ! This routine traces all extrema for a 2D (magnetic) field.
        ! The routine then traces all points in the bounded domain that satisfy:
        !   
        !       dFdx = 0
        !       dFdy = 0
        !
        ! Depending on the higher order derivatives in those points, they can then
        ! be classified as 
        !       
        !       local minimum: [d2Fdx2 d2Fdxdy; d2Fdydx d2Fdy2] is positive
        !       definite
        !       saddle point: non-positive definite
        !       local maximum:  negative definite
        !
        ! This is checked by computing the eigenvalues of the local hessian.
        ! In this case, this is relatively trivial, since the matrix is 
        ! only 2D and hence the eigenvalues can be computed analytically.  

        ! Algorithm
        !==========
        ! Actually, the algorithm is rather simple since only 2 independent
        ! variables (x, y) are present. We trace contour lines of dFdx = 0 and
        ! dFdy = 0 and then look for intersections of these lines. At
        ! intersections, we call a Newton solver to compute the exact location of
        ! dFdx = 0, dFdy = 0, and we compute locally the hessian and its
        ! eigenvalues. 

        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: when computing the optimum with the Newton solver, it is assumed
        ! that no damping strategy is necessary since the initial point should be
        ! (very) close to the optimum. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(out)      :: xe(:), ye(:), fe(:)
        integer(I8), allocatable, intent(out)   :: typee(:)
        class(ContourTracerUDT), intent(inout)  :: fieldtracer 
        type(MagneticFieldUDT), intent(in)      :: magneticField
        logical, intent(in)                     :: donewton

        ! Auxiliary
        integer(I8)                             :: ngp, nfxc, nfyc, &
            nx 
        integer(I8), allocatable, dimension(:)  :: ts1, ts2, tt
        real(R8)                                :: tempx, tempy
        real(R8), allocatable, dimension(:)     :: xg, yg, f, fx, fy, &
            tx, ty, thiseig, tf, tfxx, tfxy, tfyy
        logical                                 :: conv 
        type(RealDynamicArrayUDT)               :: xc, yc, fc 
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: fxc(:), fyc(:)
        type(PolygonSetUDT)                     :: tempps
        type(PolygonUDT), allocatable           :: fxp(:), fyp(:)
        class(ContourTracerUDT), allocatable    :: fxtracer, fytracer

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Check inputs
        if (.not. donewton) then 
            print *, 'TraceExtrema2DBox: not applying Newton solver, ' // & 
                'locations of extrema may be slightly inaccurate'
        end if

        ! Extract grid coordinates
        call fieldtracer%GetCoordinates(xg, yg)
        ngp = size(xg)

        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Trace contours
        !===============
        ! Evaluate derivatives
        allocate(f(ngp), fx(ngp), fy(ngp))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, f)
        call magneticField%interp%Evaluate(xg, yg, 1, 0, fx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, fy)

        ! Construct new tracers
        fxtracer = fieldtracer
        fytracer = fieldtracer 
        if (allocated(fxtracer%xs)) then 
            deallocate(fxtracer%xs, fxtracer%ys, fxtracer%vs, &
                fxtracer%order, fxtracer%IDs)
        end if 
        allocate(fxtracer%xs(0), fxtracer%ys(0), fxtracer%vs(0), &
            fxtracer%order(0), fxtracer%IDs(0))
        if (allocated(fytracer%xs)) then 
            deallocate(fytracer%xs, fytracer%ys, fytracer%vs, &
            fytracer%order, fytracer%IDs)
        end if 
        allocate(fytracer%xs(0), fytracer%ys(0), fytracer%vs(0), &
            fytracer%order(0), fytracer%IDs(0))
        call fxtracer%SetValues(fx)
        call fytracer%SetValues(fy)
        

        ! Trace 
        fxc = fxtracer%TraceContours([0.0_R8])
        fyc = fytracer%TraceContours([0.0_R8])

        ! Clean
        call CleanContours(fxc)
        call CleanContours(fyc)

        ! Convert all to polygons
        nfxc = size(fxc)
        nfyc = size(fyc)
        allocate(fxp(nfxc), fyp(nfyc))

        do i = 1, nfxc 
            call fxp(i)%Construct(fxc(i)%x, fxc(i)%y)
        end do 
        do i = 1, nfyc 
            call fyp(i)%Construct(fyc(i)%x, fyc(i)%y)
        end do 

        ! Write out data
        call tempps%Construct(fxp)
        call tempps%WriteData('extrema_fx_lines')
        call tempps%Construct(fyp)
        call tempps%WriteData('extrema_fy_lines')
        

        ! Compute all intersections
        !==========================
        ! Compute intersections
        !$omp parallel do default(none) schedule(dynamic) collapse(2) &
        !$omp shared(nfxc, nfyc, fxc, fyc, magneticField, xc, yc, fc, tc) &
        !$omp private(i, j, tx, ty, ts1, ts2, nx, tt, k, tempx, tempy, &
        !$omp conv, tf, tfxx, tfxy, tfyy, thiseig)
        do i = 1, nfxc
            ! Only compute intersections with other polygons
            do j = 1, nfyc
                ! Compute intersections with next polygon
                call SimplePolygonIntersections(fxc(i)%x, fxc(i)%y, &
                    fyc(j)%x, fyc(j)%y, tx, ty, ts1, ts2)
                !call PolygonIntersections(fxp(i), fyp(j), tx, ty, ts1, ts2)
                
                ! Check if found
                nx = size(tx)
            
                ! Add intersections
                allocate(tt(nx))
                do k = 1, nx 
                    ! Refine
                    call TinyNewtonSolver(tempx, tempy, conv, &
                        tx(k), ty(k), magneticField)
                     
                    ! Check 
                    if (conv) then 
                        tx(k) = tempx
                        ty(k) = tempy
                    else
                        print *, 'TraceExtrema2DBox: newton solver did ' // & 
                            'not converge, taking original estimate'
                    end if
                end do 

                ! Compute value at location and second order derivatives
                allocate(tf(nx), tfxx(nx), tfxy(nx), tfyy(nx))
                call magneticField%interp%Evaluate(tx, ty, 0, 0, tf)
                call magneticField%interp%Evaluate(tx, ty, 2, 0, tfxx)
                call magneticField%interp%Evaluate(tx, ty, 1, 1, tfxy)
                call magneticField%interp%Evaluate(tx, ty, 0, 2, tfyy)

                ! Determine order
                do k = 1, nx

                    ! Compute value and check type of extremum
                    thiseig = ComputeEigenvaluesSymmetric2by2Matrix(&
                        tfxx(k), tfyy(k), tfxy(k))
                    if (any(thiseig == 0)) then 
                        print *, 'TraceExtrema2DBox: extrema with zero ' // & 
                            'eigenvalue detected, may be improperly identified'
                    end if 
                    if (all(thiseig > 0)) then 
                        ! Local minimum
                        tt(k) = -1
                    elseif (all(thiseig < 0)) then 
                        ! Local maximum
                        tt(k) = 1
                    else
                        ! Saddle point
                        tt(k) = 0
                    end if
                end do 

                ! Append
                !$omp critical
                call xc%Append(tx)
                call yc%Append(ty)
                call fc%Append(tf)
                call tc%Append(tt)
                !$omp end critical

                ! Housekeeping
                deallocate(tt, tf, tfxx, tfyy, tfxy)
            end do 
        end do
        !$omp end parallel do 
        
        ! Extract
        !========
        xe = xc%Get()
        ye = yc%Get()
        fe = fc%Get()
        typee = tc%Get()

    end subroutine 

    ! 2D tangency points
    subroutine TraceTangencyPoints2D(xtp, ytp, typetp, Ftp, &
        boundarytracer, magneticField)

        ! Description
        !============
        ! This routine traces the tangency points on a prescribed boundary.
        ! We base ourselves solely on the
        ! discrete representation of the vessel geometry and compute the tangency
        ! points by looking at the field evaluated in the vertices of the vessel
        ! polygon(s). Here, minima and maxima are quite easily and rapidly found
        ! without having to compute any intersections between curves. To determine
        ! the type of tangency point (i.e. whether the curve bends inside or
        ! outside the domain), we check whether it's a minima or maxima along the
        ! vessel curve and whether the tangency point field value is found in the
        ! interior of the domain near the tangency point location. 

        ! The latter is a bit tricky here: since we don't explicitly know the
        ! magnetic field in the interior, we evaluate the derivative in the
        ! tangency point location and base ourselves on that. This, however, may be
        ! inaccurate since the tangency point location itself will not correspond
        ! with the location of the continuous optimum. Caution is therefore advised
        ! when the local magnetic field varies strongly locally (or when the mesh
        ! is too coarse). 


        ! Notes
        !======
        ! Note 1: the contour tracing algorithm used to trace dFdx = 0 and dFdy = 0
        ! contours is assumed to be able to deal with saddle points and should
        ! return only 'simple' polygons. Saddle points may exist even in the dFdx =
        ! 0 and dFdy = 0 fields. Therefore, we've implemented our own tracing
        ! routine in TraceContourLineStructured2D. 

        ! Note 2: it is assumed that the vessel polygon forms a closed (or a set of
        ! closed) surfaces. These polygons are oriented later on such that the face
        ! normal points inward the domain at all times. 


        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(out)      :: xtp(:), ytp(:), &
            Ftp(:)
        integer(I8), allocatable, intent(out)   :: typetp(:)
        class(ContourTracerUDT), intent(in)     :: boundarytracer
        type(MagneticFieldUDT), intent(in)      :: magneticField 

        ! Auxiliary
        integer(I8)                             :: flag
        integer(I8), allocatable                :: dvals(:), ddvals(:), &
            tv(:), extrlocind(:), tt(:)
        real(R8)                                :: fdifftol 
        real(R8), allocatable, dimension(:)     :: val, dval, tx, ty, &
            tf, nxpe, nype, nxp, nyp, normprod, dFdx, dFdy
        logical                                 :: hasbeendeleted
        logical, allocatable                    :: extrloc(:)
        type(RealDynamicArrayUDT)               :: xc, yc, fc
        type(IntegerDynamicArrayUDT)            :: tc
        type(ContourUDT), allocatable           :: bndcontours(:)
        type(PolygonUDT), allocatable           :: bndpol(:)
        type(PolygonSetUDT)                     :: bndps

        ! Loop
        integer(I8)                             :: i, j, k 

        ! Initialize
        !===========
        ! Initialize dynamic arrays
        xc = ConstructRealDynamicArray()
        yc = ConstructRealDynamicArray()
        fc = ConstructRealDynamicArray()
        tc = ConstructIntegerDynamicArray()

        ! Trace contours
        !---------------
        ! Boundary
        bndcontours = boundarytracer%TraceContours([0.0_R8])
        
        ! Sanity check
        if (size(bndcontours) == 0) then 
            ! Throw error - no boundary polygon
            call gdErrorHandler('TraceTangencyPoints2D: could not trace ' // & 
                'boundary contour, check input')
        end if 

        ! Construct boundary polygonset
        !==============================
        ! Construct polygons from boundary
        allocate(bndpol(size(bndcontours)))
        do i = 1, size(bndcontours)
            ! Construct
            call bndpol(i)%Construct(bndcontours(i)%x, bndcontours(i)%y)

            ! Check for closedness, if not -> error
            if (.not. bndpol(i)%isclosed) then 
                call gdErrorHandler('TraceTangencyPoints2D: boundary polygon ' // & 
                    'is not closed, not supported. Check input')
            end if 
        end do 

        ! Construct polygonset
        call bndps%Construct(bndpol)

        ! Orient
        call bndps%OrientNestedClosedPolygons(flag)

        ! Check if successful
        if (flag /= 0) then
            call gdErrorHandler('TraceTangencyPoints2D: could not orient ' // & 
                'boundary polygons, check input')
        end if 

        ! Compute tangency points
        !========================
        do i = 1, size(bndps%polygons)
            ! Associate for ease
            associate(p         => bndps%polygons(i))

            ! Evaluate field at vertex locations
            allocate(val(size(p%vert)))
            call magneticField%interp%Evaluate(p%x(p%vert), p%y(p%vert), &
                0, 0, val) ! assumed start and end point the same
            val = [val, val(2)] ! extend to take next edge into account

            ! Take difference
            dval = val(2:size(val)) - val(1:size(val)-1)

            ! Check where this changes sign
            allocate(dvals(size(dval)))
            where (dval > 0) dvals = 1
            where (dval <= 0) dvals = -1
            ddvals = dvals(2:size(dvals)) - dvals(1:size(dvals)-1) 

            ! Find the location and value of extrema
            extrloc = [.false., ddvals /= 0]
            allocate(tv(count(extrloc)))
            tv = pack(p%vert, extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(tv) ! should be fine
            
            ! Check if we should exclude extremum pairs based on field value
            ! difference
            fdifftol = (maxval(val) - minval(val))*tprelfieldtol
            k = 1
            hasbeendeleted = .false.
            do while (k < size(tf))
                ! Check difference
                if (abs(tf(k+1)-tf(k)) < fdifftol) then 
                    ! Remove values, such that subsequent ones can be
                    ! checked too
                    hasbeendeleted = .true.
                    extrloc(tv(k:k+1)) = .false.
                    tv = [tv(:k-1), tv(k+1:)]
                    tf = [tf(:k-1), tf(k+1:)]
                    k = 1
                else
                    k = k + 1
                end if
            end do 
            
            ! Check last 'edge'
            if (abs(tf(size(tf))-tf(1)) < fdifftol) then 
                extrloc(tv([1, size(tf)])) = .false.
                tv = pack(p%vert, extrloc)
                tf = val(tv)
                hasbeendeleted = .true.
            end if
            
            ! Issue message
            if (hasbeendeleted) then 
                print *, 'TraceTangencyPoints2DBoxDiscrete: some ' // & 
                    'tangency points were deleted based on their field ' // & 
                    'values as they are very close together'
            end if

            ! Recompute points
            tv = pack(p%vert, extrloc)
            extrlocind = tv 
            extrlocind = pack([(k, k = 1, size(extrloc))], extrloc)
            tx = p%x(tv)
            ty = p%y(tv)
            tf = val(tv) ! should be fine

            ! Compute product between normal at vertex and magnetic field
            allocate(dFdx(size(tf)), dFdy(size(tf)))
            call magneticField%interp%Evaluate(tx, ty, 1, 0, dFdx)
            call magneticField%interp%Evaluate(tx, ty, 0, 1, dFdy)

            nxpe = [p%nx, p%nx(1)]
            nype = [p%ny, p%ny(1)]
            nxp = 0.5*(nxpe(extrlocind-1) + nxpe(extrlocind))
            nyp = 0.5*(nype(extrlocind-1) + nype(extrlocind))
            normprod = dFdx*nxp + dFdy*nyp

            ! Determine type and add
            call xc%Append(tx)
            call yc%Append(ty)
            call fc%Append(tf)
            allocate(tt(size(tf)))
            tt = 0
            do j = 1, size(extrlocind)

                ! Determine type
                if ((val(extrlocind(j)-1) < tf(j)) .and. &
                    (val(extrlocind(j)+1) < tf(j))) then 

                    ! Local maximum
                    if (normprod(j) < 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = 1
                    else
                        ! Curve bends inwards
                        tt(j) = -1
                    end if
                else
                    ! Local minimum
                    if (normprod(j) > 0) then 
                        ! Curve bends outwards of the domain
                        tt(j) = 1
                    else
                        ! Curve bends inwards
                        tt(j) = -1
                    end if 
                end if
            end do 

            ! Append
            call tc%Append(tt)

            ! Housekeeping
            deallocate(val, dvals, dFdx, dFdy, tt, tv)
            end associate

        end do 

        ! Extract values
        !===============
        xtp = xc%Get()
        ytp = yc%Get()
        Ftp = fc%Get()
        typetp = tc%Get()



    end subroutine

    !------------------------------------------------------------------!
    !                    TOPOLOGICAL MESH OPERATORS                    !
    !------------------------------------------------------------------!

    ! Initializers
    subroutine InitializeTopologicalMesh(topomesh)

        ! Description
        !============
        ! Initialize the topological mesh substructures

        ! Declare variables
        !==================
        class(TopomeshUDT)      :: topomesh 

        ! Initialize
        !===========
        ! number of distinct flux surfaces
        topomesh%nfs = 0
        topomesh%fsfval = ConstructRealDynamicArray()
        topomesh%fsID   = ConstructIntegerDynamicArray()

        ! Substructures (tubes are initialized empty)
        call topomesh%vert%Initialize()
        call topomesh%face%Initialize()
        call topomesh%cell%Initialize()
        call topomesh%tube%Initialize(0, 0, 0, 0, 0, 0, 0)

    end subroutine

    subroutine InitializeTopologicalMeshVertex(tpvert, nv)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshVertUDT)      :: tpvert 
        integer(I8), intent(in), optional   :: nv 

        ! Initialize
        !===========
        if (present(nv)) then 
            tpvert%ntot = nv 
        else
            tpvert%ntot = 0
        end if 
        if (allocated(tpvert%ID)) then 
            ! Assume all allocated
            deallocate(tpvert%ID, tpvert%x, tpvert%y, tpvert%fval, &
                tpvert%type, tpvert%fsID, tpvert%BV)
        end if 
        associate(ntot  => tpvert%ntot)
        allocate(tpvert%ID(ntot), tpvert%x(ntot), tpvert%y(ntot), &
            tpvert%fval(ntot), tpvert%type(ntot), tpvert%fsID(ntot), &
            tpvert%BV(ntot))

        end associate

    end subroutine 

    subroutine InitializeTopologicalMeshFace(tpface, nf, nfc)

        ! Description
        !============
        ! Initialize the topomesh vertex structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshFaceUDT)      :: tpface
        integer(I8), intent(in), optional   :: nf, nfc

        ! Initialize
        !===========
        if (present(nf)) then 
            tpface%ntot = nf
        else
            tpface%ntot = 0
        end if
        if (present(nfc)) then 
            tpface%nfc = nfc
        else 
            tpface%nfc = 0
        end if 
        associate(ntot  => tpface%ntot, ntotc   => tpface%nfc)
        if (allocated(tpface%ID)) then 
            ! assume all allocated
            deallocate(tpface%ID, tpface%vert, tpface%cell, &
                tpface%fsID, tpface%x, tpface%y, &
                tpface%pol, tpface%type, tpface%cellP, tpface%label, &
                tpface%BF)
        end if 
        allocate(tpface%ID(ntot), tpface%vert(ntot, 2), tpface%cell(ntotc), &
            tpface%fsID(ntot), tpface%x(ntot), tpface%y(ntot), &
            tpface%pol(ntot), tpface%type(ntot), tpface%cellP(ntot, 2), &
            tpface%label(ntot), tpface%BF(ntot))

        ! Housekeeping
        end associate

    end subroutine

    subroutine InitializeTopologicalMeshCell(tpcell, nc, ncv, ncf)

        ! Description
        !============
        ! Initialize the topomesh cell structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshCellUDT)      :: tpcell
        integer(I8), intent(in), optional   :: nc, ncv, ncf

        ! Initialize
        !===========
        if (present(nc)) then 
            tpcell%ntot = nc
        else
            tpcell%ntot = 0
        end if 
        if (present(ncv)) then 
            tpcell%ntotv = ncv 
        else 
            tpcell%ntotv = 0
        end if 
        if (present(ncf)) then 
            tpcell%ntotf = ncf 
        else
            tpcell%ntotf = 0
        end if 

        if (allocated(tpcell%ID)) then 
            ! Assume all allocated
            deallocate(tpcell%ID, tpcell%vert, tpcell%vertP, tpcell%face,  &
                tpcell%faceP, tpcell%flags)
        end if 
        associate(ntot  => tpcell%ntot, ntotf   => tpcell%ntotf, ntotv  => tpcell%ntotv)
        allocate(tpcell%ID(ntot), tpcell%vert(ntotv), tpcell%vertP(ntot, 2), &
            tpcell%face(ntotf),  tpcell%faceP(ntot, 2), tpcell%flags(ntot))
        end associate

    end subroutine

    subroutine InitializeTopologicalMeshTube(tptube, ntot, nface, ncell, &
        nf1, nf2, nv1, nv2)

        ! Description
        !============
        ! Initialize the topomesh tube structure (simply empty arrays)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshTubeUDT)      :: tptube 
        integer(I8), intent(in)     :: ntot, nface, ncell, nv1, nv2, &
            nf1, nf2

        ! Initialize
        !===========
        tptube%ntot = ntot
        tptube%nface = nface
        tptube%ncell = ncell
        if (allocated(tptube%cell)) then 
            ! assume all allocated
            deallocate(tptube%cell, tptube%face, &
                tptube%cellP, tptube%faceP, tptube%isclosed, tptube%bndf1, tptube%bndf2, &
                tptube%bndv1, tptube%bndv2, tptube%bndf1P, &
                tptube%bndf2P, tptube%bndv1P, tptube%bndv2P, tptube%hfside)
        end if 

        allocate(tptube%cell(ncell), tptube%face(nface), &
            tptube%cellP(ntot, 2), tptube%faceP(ntot, 2), &
            tptube%isclosed(ntot), tptube%bndf1(nf1), tptube%bndf2(nf2), &
            tptube%bndv1(nv1), tptube%bndv2(nv2), tptube%bndf1P(ntot, 2), &
            tptube%bndf2P(ntot, 2), tptube%bndv1P(ntot, 2), tptube%bndv2P(ntot, 2), &
            tptube%hfside(ntot))

    end subroutine

    ! Deallocators
    subroutine DeallocateTopologicalMeshTube(tptube)

        ! Description
        !============
        ! Soft deallocation

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshTubeUDT)          :: tptube 

        ! Deallocate
        !===========
        tptube%ntot = 0
        tptube%nface = 0
        tptube%ncell = 0
        if (allocated(tptube%face)) then 
            deallocate(tptube%face)
        end if 
        if (allocated(tptube%faceP)) then 
            deallocate(tptube%faceP)
        end if 
        if (allocated(tptube%cell)) then 
            deallocate(tptube%cell)
        end if 
        if (allocated(tptube%cellP)) then 
            deallocate(tptube%cellP)
        end if 
        if (allocated(tptube%isclosed)) then 
            deallocate(tptube%isclosed)
        end if 
        if (allocated(tptube%bndf1)) then 
            deallocate(tptube%bndf1)
        end if 
        if (allocated(tptube%bndf2)) then 
            deallocate(tptube%bndf2)
        end if 
        if (allocated(tptube%bndv1)) then 
            deallocate(tptube%bndv1)
        end if 
        if (allocated(tptube%bndv2)) then 
            deallocate(tptube%bndv2)
        end if 
        if (allocated(tptube%bndf1P)) then 
            deallocate(tptube%bndf1P)
        end if 
        if (allocated(tptube%bndf2P)) then 
            deallocate(tptube%bndf2P)
        end if 
        if (allocated(tptube%bndv1P)) then 
            deallocate(tptube%bndv1P)
        end if 
        if (allocated(tptube%bndv2P)) then 
            deallocate(tptube%bndv2P)
        end if 
        if (allocated(tptube%hfside)) then 
            deallocate(tptube%hfside)
        end if  

    end subroutine

    ! Vertex addition
    subroutine AddTopologicalMeshVertex(topomesh, x, y, F, t, fsID)

        ! Description
        !============
        ! Add a vertex to the topological mesh. Not optimized for 
        ! memory usage (e.g. less allocate/deallocate), but 
        ! shouldn't be an issue since topomeshes are usually small. Note
        ! that no interconnection data is updated!

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        real(R8), intent(in)                :: x, y, F
        integer(I8), intent(in)             :: t, fsID

        ! Auxiliary
        real(R8)                            :: tF 

        ! Check 
        !======
        ! If the flux surface ID was already added, overwrite the flux
        ! value - this should be the same (but discretely this often
        ! doesn't hold) 
        tF = F 
        if (fsID /= 0) then 
            if ((topomesh%fsID%Size() > 0)) then 
                if (any(topomesh%fsID%Get() == fsID)) then 
                    ! Already exists, overwrite value
                    tF = topomesh%fsfval%Get(fsID)
                else
                    ! Does not yet exist, add 
                    call topomesh%fsID%Append(fsID)
                    call topomesh%fsfval%Set(fsID, F)
                end if
            else
                ! Does not yet exist, add 
                call topomesh%fsID%Append(fsID)
                call topomesh%fsfval%Set(fsID, F)
            end if 
        end if 


        ! Concatenate
        !============
        topomesh%vert%ntot = topomesh%vert%ntot + 1
        topomesh%vert%x = [topomesh%vert%x, x]
        topomesh%vert%y = [topomesh%vert%y, y]
        topomesh%vert%fval = [topomesh%vert%fval, tF]
        topomesh%vert%type = [topomesh%vert%type, t]
        topomesh%vert%ID = [topomesh%vert%ID, topomesh%vert%ntot]
        topomesh%vert%fsID = [topomesh%vert%fsID, fsID]

    end subroutine

    ! Vertex face addition
    subroutine AddTopologicalMeshVertexFaces(topomesh)

        ! Description
        !============
        ! This function adds the faces of each vertex and sorts them in a uniquely
        ! defined direction which is the same for all faces. It is assumed that all
        ! faces and vertices of the topology are adequately defined and that each
        ! face has a starting and end point etc (so basically after all
        ! intersections and topology faces have been added). 

        ! Algorithm
        !==========
        ! 1) For each vertex, find all faces that have this vertex and store them
        ! in vert.face and vert.faceP (list and pointer)
        ! 2) For each vertex, sort this list by doing the following steps:
        !   2.1) For each face, get the closest point of that face to the current
        !   vertex (but with distance > 1e-12 to hedge for numerical bullshit)
        !   2.2) Take one face as reference, compute the angle of all other faces
        !   w.r.t. that first face. 
        !   2.3) Sort the faces in sequence of increasing angle
        !   2.4) Check the cross product of each pair of consecutive faces. All
        !   cross-products should have the same sign. If it is positive, keep it
        !   like that. If it is negative, reverse the sorting direction. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: nfpv, tf, sortind, &
            tfsorted , fID
        real(R8)                                :: vx, vy
        real(R8), allocatable, dimension(:)     :: npx, npy, dx, dy,&
            theta
        logical, allocatable                    :: test(:)

        ! Loop
        integer(I8)                             :: i, j, k

        ! Initialize
        !===========
        ! Do some checks
        if (any(topomesh%face%vert(:, 1) == 0) .or. any(topomesh%face%vert(:, 2) == 0)) then  
            call gdErrorHandler('AddTopologicalMeshVertexFaces: faces without vertex ' // & 
                'encountered. First, clean up the topological mesh before calling this routine')
        end if 

        ! Get faces of each vertex
        !-========================
        ! Compute for each vertex how many faces they would have
        allocate(nfpv(topomesh%vert%ntot))
        nfpv = 0
        do i = 1, topomesh%face%ntot
            nfpv(topomesh%face%vert(i, :)) = nfpv(topomesh%face%vert(i, :)) + 1
        end do

        ! Initialize 
        if (allocated(topomesh%vert%face)) then 
            deallocate(topomesh%vert%face)
        end if 
        if (allocated(topomesh%vert%faceP)) then 
            deallocate(topomesh%vert%faceP)
        end if 
        allocate(topomesh%vert%face(sum(nfpv)))

        ! Construct the face pointer
        allocate(topomesh%vert%faceP(topomesh%vert%ntot, 2))
        topomesh%vert%faceP(:, 2) = nfpv 
        topomesh%vert%faceP(1, 1) = 1
        do i = 2, topomesh%vert%ntot 
            topomesh%vert%faceP(i, 1) = topomesh%vert%faceP(i-1, 1) + nfpv(i-1)
        end do 

        ! Get face list (unsorted)
        allocate(fID(topomesh%face%ntot))
        fID = [(k, k = 1, topomesh%face%ntot)]
        do i = 1, topomesh%vert%ntot
            ! Get face indices
            test = (topomesh%face%vert(:, 1) == i) .or. (topomesh%face%vert(:, 2) == i) 
            tf = pack(fID, test)
            
            ! Add
            topomesh%vert%face(topomesh%vert%faceP(i, 1):topomesh%vert%faceP(i, 1)+topomesh%vert%faceP(i, 2)-1) = tf
        end do

        ! Sort faces
        !===========
        ! Sort
        do i = 1, topomesh%vert%ntot
            ! Unpack this vertex
            vx = topomesh%vert%x(i)
            vy = topomesh%vert%y(i)

            ! Get faces
            tf = GetTMVertFace(topomesh%vert, i)
            
            ! Get the next point of each face 
            allocate(npx(size(tf)), npy(size(tf)))
            do j = 1, size(tf)
                if (topomesh%face%vert(tf(j), 1) == i) then  ! first vertex
                    npx(j) = topomesh%face%x(tf(j))%Get(2)
                    npy(j) = topomesh%face%y(tf(j))%Get(2)
                elseif (topomesh%face%vert(tf(j), 2) == i) then  ! second vertex
                    npx(j) = topomesh%face%x(tf(j))%Get(topomesh%face%x(tf(j))%size()-1)
                    npy(j) = topomesh%face%y(tf(j))%Get(topomesh%face%y(tf(j))%size()-1)
                else
                    call gdErrorHandler('AddTopologicalMeshVertexFaces: ' // & 
                        'This should not be happening and is a bug!')
                end if 
            end do 
            
            ! Compute angles (w.r.t. horizontal axis)
            dx = npx - vx
            dy = npy - vy
            theta = atan2(dy, dx)
            
            ! Sort
            allocate(sortind(size(theta)))
            call Sort(theta, ind=sortind)
            tfsorted = tf(sortind)
            deallocate(sortind)
            
            ! Add
            topomesh%vert%face(topomesh%vert%faceP(i, 1):&
                topomesh%vert%faceP(i, 1)+topomesh%vert%faceP(i, 2)-1) = tfsorted

            ! Housekeeping
            deallocate(npx, npy)
            
        end do


    end subroutine

    ! Face addition
    subroutine AddTopologicalMeshFace(topomesh, facevert, x, y, &
        t, fsID, F)

        ! Description
        !============
        ! Add a face to the topological mesh. Not optimized for 
        ! memory usage (e.g. less allocate/deallocate), but 
        ! shouldn't be an issue since topomeshes are usually small. Note
        ! that no interconnection data is updated!

        ! Note: to hedge for any issues in other routines, it is checked
        ! whether subsequent coordinates are coinciding to disttol 
        ! precision. These are removed from the x, y coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        type(RealDynamicArrayUDT)           :: x, y 
        integer(I8), intent(in)             :: t, facevert(1:2), fsID
        real(R8), intent(in)                :: F 

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :)
        real(R8)                            :: tF 
        real(R8), allocatable               :: xt(:), yt(:), dxt(:), &
            dyt(:)
        logical, allocatable                :: delind(:)
        type(PolygonUDT)                    :: tp 

        ! Check if flux surface exists
        !=============================
        tF = F 
        if (fsID /= 0) then 
            if ((topomesh%fsID%Size() > 0)) then 
                if (any(topomesh%fsID%Get() == fsID)) then 
                    ! Already exists, overwrite value
                    tF = topomesh%fsfval%Get(fsID)
                else
                    ! Does not yet exist, add 
                    call topomesh%fsID%Append(fsID)
                    call topomesh%fsfval%Set(fsID, F)
                end if
            else
                ! Does not yet exist, add 
                call topomesh%fsID%Append(fsID)
                call topomesh%fsfval%Set(fsID, F)
            end if 
        end if 

        ! Hedge for duplicate points
        !===========================
        ! Detect duplicate points
        xt = x%Get()
        yt = y%Get()
        dxt = xt(2:size(xt)) - xt(1:size(xt)-1)
        dyt = yt(2:size(yt)) - yt(1:size(yt)-1)
        allocate(delind(size(xt)))
        delind = .false. 
        delind(2:size(delind)) = ((abs(dxt) <= disttol) .and. (abs(dyt)) <= disttol)
        if (any(delind)) then 
            x = ConstructRealDynamicArray(pack(xt, .not. delind))
            y = ConstructRealDynamicArray(pack(yt, .not. delind))
        end if 

        ! Make sure that final and start point are equal to vertex 
        ! coordinates if facevert is not zero
        if (x%size() < 2) then 
            print *, 'size is smaller than 2'
        end if 
        if (facevert(1) /= 0) then 
            call x%Set(1, topomesh%vert%x(facevert(1)))
            call y%Set(1, topomesh%vert%y(facevert(1)))
        end if 
        if (facevert(2) /= 0) then 
            call x%Set(x%Size(), topomesh%vert%x(facevert(2)))
            call y%Set(y%Size(), topomesh%vert%y(facevert(2)))
        end if 

        !  Construct new polygon
        call tp%Construct(x%Get(), y%Get())

        ! Ensure minimal amount of vertices
        if (tp%nv < 4) then ! 4 should be adequate to represent also closed polygons
            call tp%Refine(0.0_R8, 4)
            x = ConstructRealDynamicArray(tp%x(tp%vert))
            y = ConstructRealDynamicArray(tp%y(tp%vert))
        end if 


        ! Concatenate
        !============
        topomesh%face%ntot = topomesh%face%ntot + 1
        tmp = topomesh%face%vert 
        deallocate(topomesh%face%vert)
        allocate(topomesh%face%vert(topomesh%face%ntot, 2))
        topomesh%face%vert(1:topomesh%face%ntot-1, :) = tmp 
        topomesh%face%vert(topomesh%face%ntot, :) = facevert 
        topomesh%face%x = [topomesh%face%x, x]
        topomesh%face%y = [topomesh%face%y, y]
        topomesh%face%fsID = [topomesh%face%fsID, fsID]
        topomesh%face%type = [topomesh%face%type, t]
        topomesh%face%ID = [topomesh%face%ID, topomesh%face%ntot]
        topomesh%face%pol = [topomesh%face%pol, tp]

    end subroutine

    ! Face cell addition
    subroutine AddTopologicalMeshFaceCells(topomesh)

        ! Description
        !============
        ! Add the cell neighbours of a face using the pointer way of working. It is
        ! assumed that all cells etc have been constructed correctly and that other
        ! basic topology information is available. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh   

        ! Auxiliary
        integer(I8)                             :: ind
        integer(I8), allocatable, dimension(:)  :: ncpf, tf, fc

        ! Loop
        integer(I8)                             :: i, j

        ! Add faces
        !==========
        ! Initialize
        if (allocated(topomesh%face%cellP)) then 
            deallocate(topomesh%face%cellP)
        end if 
        if (allocated(topomesh%face%cell)) then 
            deallocate(topomesh%face%cell)
        end if 
        allocate(topomesh%face%cellP(topomesh%face%ntot, 2))

        ! Compute number of cells per face
        allocate(ncpf(topomesh%face%ntot))
        ncpf = 0
        do i = 1, topomesh%cell%ntot
            ! Get cell faces
            tf = GetTMCellFace(topomesh%cell, i)
            ncpf(tf) = ncpf(tf) + 1
        end do

        ! Set cellP
        topomesh%face%cellP(:, 2) = ncpf
        topomesh%face%cellP(1, 1) = 1
        do i = 2, topomesh%face%ntot 
            topomesh%face%cellP(i, 1) = topomesh%face%cellP(i-1, 1) + & 
                topomesh%face%cellP(i-1, 2)
        end do 

        ! Set cells of each face
        allocate(topomesh%face%cell(sum(ncpf)), fc(topomesh%face%ntot))
        fc = 0
        topomesh%face%cell = 0
        do i = 1, topomesh%cell%ntot
            ! Get cell faces
            tf = GetTMCellFace(topomesh%cell, i);
            do j = 1, size(tf)
                ! Add
                ind = topomesh%face%cellP(tf(j), 1) + fc(tf(j))
                topomesh%face%cell(ind) = i
                
                ! Update counter
                fc(tf(j)) = fc(tf(j)) + 1
            end do 
        end do 

    end subroutine 

    ! Cell addition 
    subroutine AddTopologicalMeshCells(topomesh)

        ! Description
        !============
        ! This routine forms cells of the given topology mesh with consistent data
        ! of faces and vertices. It is assumed that all faces are 'simple', i.e.
        ! they do not start and end in the same point (in that case, points
        ! should've been introduced before). The cells can be of arbitrary topology
        ! (i.e. triangles, quads, ...) - if this is not desired, it should be
        ! checked afterwards.

        ! The following cell types are allowed: 
        ! - regular cells: cells with all unique boundaries that do not appear
        ! twice for this cell 
        ! - Disc-type cells: cells that have a single boundary appearing twice, but
        ! where this single boundary has one vertex that only has one face.

        ! All other cell types will not be properly identified by this routine and
        ! will very likely lead to undesired output or even crashing. If necessary,
        ! some other cell types may be supported in the future by making the stop
        ! condition for a cell more stringent (i.e. not just stopping if we find
        ! the same face, but if we also arrive at the face in the same way as we
        ! started - this would allow cells with duplicate boundaries as well). 

        ! Algorithm
        !==========
        ! We assume that the following data is readily available:
        ! - face indices of vertices, sorted either clockwise or counter-clockwise,
        ! but consistently the same for each vertex. This allows us to uniquely
        ! define the left and right neighbor face when determining the next one.
        ! This should be computed beforehand with the
        ! 'AddTopologicalMeshVertexFaces' routine.
        ! - Face identifiers that are accurate and up to date: at least boundary
        ! faces must be indicated to be type 3. No closed faces are assumed to
        ! exist. 

        ! Furthermore, we assume that the ensemble of vertices and faces results in
        ! a non-overlapping partition of the bounded 2D domain, and that the union
        ! of all faces with type 3 represents a closed boundary polygon (or
        ! multiple closed polygons) that bound the domain.

        ! Declare variables
        !==================
        ! Arguments 
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary 
        integer(I8)                             :: tf, startface, & 
            turndirection, tfv(1:2), tvind, tv, startvert, &
            starttvind, nf, nfv(1:2), fu
        integer(I8), allocatable, dimension(:)  :: fc, disccellvert, &
            nfvfn, tcf, tcv, faceneig1, faceneig2, faceneig

        logical                                 :: istfv(1:2)
        logical, allocatable                    :: hasturned1(:, :), &
            hasturned2(:, :), donotstartfromface(:)

        type(IntegerDynamicArrayUDT), allocatable   :: cellvert(:), &
            cellface(:)
        type(IntegerDynamicArrayUDT)                :: thiscellvert, &
            thiscellface

        ! Loop 
        integer(I8)                             :: i, cc 

        ! Initialize
        !===========
        ! Unpack
        associate( & 
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            cell        => topomesh%cell    &
            )

        ! Allocate
        allocate(fc(face%ntot), hasturned1(face%ntot, 2), &
            hasturned2(face%ntot, 2), donotstartfromface(face%ntot), &
            cellvert(0), cellface(0))

        ! Initialize face counters
        where (face%BF) 
            fc = 1 ! boundary faces: only once
        elsewhere 
            fc = 2
        end where 
            
        ! Initialize turn checkers
        hasturned1 = .false.
        hasturned2 = .false.

        ! Initialize cell counter
        cc = 0

        ! Check if there are any disc-type cells
        disccellvert = findloc(vert%faceP(:, 2), 1_I8)

        ! Check if there are any faces with only one adjacent face on each side -
        ! these faces shouldn't be started from, as one cannot determine the
        ! turning direction (and several other faces should remain that can be
        ! started from)
        donotstartfromface = .false. 
        do i = 1, face%ntot
            nfvfn = vert%faceP(face%vert(i, :), 2)
            if (all(nfvfn == 2)) then 
                donotstartfromface(i) = .true.
            end if 
        end do 

        ! Loop
        !=====
        do while (.true.)
            
            ! Initialize cell faces & vertices
            allocate(tcv(0), tcf(0))

            ! Find the next face (any next internal face)
            tf = findloc( (fc > 0) .and. (.not. face%BF) .and. &
                (.not. donotstartfromface), .true., 1)
            
            ! Check
            if (tf == 0) then 
                ! Check
                if (any(fc > 0)) then
                    ! Check if there are any internal faces to begin with
                    if (count(.not. face%BF) == 0) then 
                        print *, 'AddTopologicalMeshCells: case without any ' // & 
                            'internal faces detected'

                        ! Find the next face (any next face)
                        tf = findloc( (fc > 0)  , .true., 1)
                    else
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'could not find next face')
                    end if 
                else
                    ! All faces added, exit
                    exit 
                end if 
            end if 
            
            ! Do not subtract a counter - we need to end up in this face again.
            ! Also, don't add, we do this later on
            
            ! Set starting face for this cell
            startface = tf
            
            ! Turning direction (0: not found, 1: first neighbour, 2: second
            ! neighbour)
            turndirection = 0
            
            ! Get neighbouring faces in correct order
            tfv = face%vert(tf, :)
            faceneig1 = GetTMVertFaceNeig(vert, tfv(1), tf);
            faceneig2 = GetTMVertFaceNeig(vert, tfv(2), tf);
            
            ! Sanity checks
            if ((size(faceneig1) == 0) .or. (size(faceneig2) == 0)) then 
                ! No neighbours found
                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                    'could not find neighbouring faces, something wrong ' //& 
                    'with topological mesh construction. Check input')
            end if 
            if (all(faceneig1 == tf) .and. all(faceneig2 == tf)) then 
                ! Isolated face
                call gdErrorHandler('AddTopologicalMeshCells: isolated ' // & 
                    'face found, check input')
            end if
            
            ! Set the starting vertex index (if boundary vertices, needs checks)
            tvind = 0
            if ((.not. all(faceneig1 == tf)) .and. (.not. vert%BV(tfv(1)))) then 
                tvind = 1
            elseif ((.not. all(faceneig2 == tf)) .and. (.not. vert%BV(tfv(2)))) then 
                tvind = 2
            end if 

            ! If none found,  check if the first vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig1) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig1(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig1(1) == faceneig1(2)) then 
                    tvind = 1
                else
                    ! Check if any face neighbour can be taken 
                    if (any((fc(faceneig1) > 0))) then 
                        tvind = 1
                    end if 
                end if 
            end if 

            ! If none found, check if the second vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig2) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig2(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig2(1) == faceneig2(2)) then 
                    tvind = 2
                else
                    ! Check if any face neighbour can be taken
                    if (any((fc(faceneig2) > 0))) then  
                        tvind = 2
                    end if 
                end if 
            end if 
            
            ! Sanity check
            if (tvind == 0) then 
                ! No starting vertex found - possibly dangling face?
                call gdErrorHandler('AddTopologicalMeshCells: could not ' // & 
                    'find next cell although initial face was found - ' // & 
                    'possible dangling face detected, check input')
            end if
            
            ! Set current vertex
            tv = face%vert(tf, tvind)
            startvert = tv
            starttvind = tvind
        
            ! Loop 
            !=====
            do while (.true.)
                
                ! Get neighbouring faces in correct order
                faceneig = GetTMVertFaceNeig(vert, tv, tf)
                        
                ! Sanity checks
                if (size(faceneig) == 0) then 
                    ! No neighbours found
                    call gdErrorHandler('AddTopologicalMeshCells: could ' // & 
                        'not find neighbouring faces, something wrong ' // & 
                        'with topological mesh construction. Check input')
                end if 
                if (all(fc(faceneig) <= 0)) then 
                    ! No neighbours with counter left
                    print *, 'vertex: ', tv
                    call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                    call gdErrorHandler('AddTopologicalMeshCells: all ' // & 
                        'neighbouring faces cannot be taken anymore, ' // &
                        'faces do not seem to form cell')
                end if 
                
                ! Find the next face
                if (turndirection /= 0) then 
                    ! We have a turn direction, so we can only check if we should
                    ! throw errors
                    nf = faceneig(turndirection)
                    
                    ! Check counter
                    if (fc(nf) <= 0) then 
                        print *, 'vertex: ', tv
                        call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                        open (action='write', file='./output/cell_error.dat', newunit=fu, &
                            status='unknown')
                        write (fu, *) 'cells that have been constructed. ' // & 
                            'per cell, each row has: vertices, face, face type'
                        do i = 1, size(cellvert)
                            write(fu, *) 'cell: ', i
                            write(fu, *) cellvert(i)%Get()
                            write(fu, *) cellface(i)%Get()
                            write(fu, *) topomesh%face%type(cellface(i)%Get())
                        end do 
                        close(fu)
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'next face is forced by turning direction ' // & 
                            'but is not available')
                    end if 
                    
                    ! If we passed this, we should check the vertices
                    nfv = face%vert(nf, :)
                    istfv = nfv == tv
                    if (.not.any(istfv)) then 
                        ! Current vertex is not found in the next face, this should
                        ! not be possible
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'next face does not have current vertex, ' // & 
                            'check input')
                    end if 
                    if (all(istfv)) then 
                        ! Next face is a face that starts and ends in the same
                        ! vertex - not supported
                        call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                            'face detected with same start and end vertex, ' // & 
                            'not supported')
                    end if 
                    
                    ! Add and update
                    tcf = [tcf, nf]
                    tcv = [tcv, tv]
                    tf = nf
                    if (.not. istfv(1)) then 
                        tv = nfv(1)
                    else
                        tv = nfv(2)
                    end if 
                    
                    ! Update counter
                    fc(tf) = fc(tf) - 1
                    
                    ! Set that this direction can't be turned in anymore
                    ! from neither side for this face
                    if ((startface == nf) .and. (tv == startvert)) then 
                        ! Just break, turning direction etc already adjusted
                        ! before.
                        exit
                    else
                        if (turndirection == 1) then 
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            else 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            end if 
                        elseif (turndirection == 2) then 
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            else 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            end if 
                        else
                            call gdErrorHandler('AddTopologicalMEshCells: ' // & 
                                'bug detecetd when adjusting turning direction')
                        end if 
                    end if 
                    
                else
                    ! We don't have a turning direction yet. Check the current
                    ! neighbours
                    if (faceneig(1) == faceneig(2)) then 
                        ! We have to take this face if we can, or we should error.
                        ! Turning direction cannot be determined
                        if (fc(faceneig(1)) <= 0) then 
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'next face is only possible face but ' // & 
                                'cannot be taken due to counter being zero')
                        end if             
                        
                        ! Take face
                        nf = faceneig(1)
                        
                        ! If we passed this, we should check the vertices
                        nfv = face%vert(nf, :)
                        istfv = nfv == tv
                        if (.not. any(istfv)) then 
                            ! Current vertex is not found in the next face, this should
                            ! not be possible
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'next face does not have current vertex, ' // & 
                                'check input')
                        end if 
                        if (all(istfv)) then 
                            ! Next face is a face that starts and ends in the same
                            ! vertex - not supported
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'face detected with same start and end ' // & 
                                'vertex, not supported')
                        end if

                        ! Add and update
                        tcf = [tcf, nf]
                        tcv = [tcv, tv]
                        tf = nf
                        if (.not. istfv(1)) then 
                            tv = nfv(1)
                        else
                            tv = nfv(2)
                        end if 
                        
                        ! Update counter
                        fc(tf) = fc(tf) - 1
                        
                        ! Is the next face the start face? If so, exit
                        if ((startface == nf) .and. (tv == startvert)) then 
                            exit 
                        end if 
                    else
                        ! Check which one has a positive counter (at least one,
                        ! cause we already passed a check for this)
                        if ((fc(faceneig(1)) > 0) .and. .not. hasturned1(startface, starttvind)) then 
                            ! Take face
                            nf = faceneig(1)
                            
                            ! If we passed this, we should check the vertices
                            nfv = face%vert(nf, :)
                            istfv = nfv == tv
                            if (.not. any(istfv)) then 
                                ! Current vertex is not found in the next face, this should
                                ! not be possible
                                call gdErrorHandler('AddTopologicalMeshCells:  ' // & 
                                    'next face does not have current vertex,  ' // & 
                                    'check input')
                            end if
                            if (all(istfv)) then 
                                ! Next face is a face that starts and ends in the same
                                ! vertex - not supported
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'face detected with same start and ' // & 
                                    'end vertex, not supported')
                            end if 
                            
                            ! Add and update
                            tcf = [tcf, nf]
                            tcv = [tcv, tv]
                            tf = nf
                            if (.not. istfv(1)) then 
                                tv = nfv(1)
                            else
                                tv = nfv(2)
                            end if 
                            
                            ! Update counter
                            fc(tf) = fc(tf) - 1
                            
                            ! Set turn direction
                            turndirection = 1
                            
                            ! Set that this direction can't be turned in anymore
                            ! from neither side
                            hasturned1(startface, starttvind) = .true.
                            if (starttvind == 1) then 
                                hasturned2(startface, 2) = .true.
                            else
                                hasturned2(startface, 1) = .true.
                            end if
                            
                            ! Remove the turn direction from the next face (otherwise this is skipped...)
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            else 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            end if 
                            
                            ! Is the next face the start face? If so, exit
                            if ((startface == nf) .and. (tv == startvert)) then 
                                exit 
                            end if 
                            
                        elseif ((fc(faceneig(2)) > 0) .and. .not. hasturned2(startface, starttvind)) then 
                            
                            ! Take face
                            nf = faceneig(2)
                            
                            ! If we passed this, we should check the vertices
                            nfv = face%vert(nf, :)
                            istfv = nfv == tv;
                            if (.not. any(istfv)) then 
                                ! Current vertex is not found in the next face, this should
                                ! not be possible
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'next face does not have current vertex, ' // & 
                                    'check input')
                            end if
                            if (all(istfv)) then 
                                ! Next face is a face that starts and ends in the same
                                ! vertex - not supported
                                call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                    'face detected with same start and ' // & 
                                    'end vertex, not supported')
                            end if
                            
                            ! Add and update
                            tcf = [tcf, nf]
                            tcv = [tcv, tv]
                            tf = nf
                            if (.not. istfv(1)) then 
                                tv = nfv(1)
                            else
                                tv = nfv(2)
                            end if 
                            
                            ! Update counter
                            fc(tf) = fc(tf) - 1
                            
                            ! Set turn direction
                            turndirection = 2
                            
                            ! Set that this direction can't be turned in anymore
                            ! from neither side
                            hasturned2(startface, starttvind) = .true.
                            if (starttvind == 1) then 
                                hasturned1(startface, 2) = .true.
                            else
                                hasturned1(startface, 1) = .true.
                            end if 

                            ! Remove the turn direction from the next face (otherwise this is skipped...)
                            if (face%vert(tf, 1) == tv) then 
                                hasturned1(tf, 2) = .true. 
                                hasturned2(tf, 1) = .true. 
                            else 
                                hasturned1(tf, 1) = .true. 
                                hasturned2(tf, 2) = .true. 
                            end if 
                            
                            ! Is the next face the start face? If so, exit
                            if ((startface == nf) .and. (tv == startvert)) then 
                                exit 
                            end if 
                            
                        else
                            call gdErrorHandler('AddTopologicalMeshCells: ' // & 
                                'could not find next face since no ' // &
                                'combination of available start direction ' // &
                                'and free face found')
                        end if 
                    end if 
                end if 
            end do 
            
            ! Add found cell to the structure
            cc = cc + 1
            thiscellvert = ConstructIntegerDynamicArray(tcv)
            thiscellface = ConstructIntegerDynamicArray(tcf)
            cellvert = [cellvert, thiscellvert]
            cellface = [cellface, thiscellface]

            ! Housekeeping
            deallocate(tcv, tcf)
        end do 

        ! Add to the topology mesh
        if (allocated(topomesh%cell%vert)) then 
            deallocate(topomesh%cell%vert)
        end if 
        if (allocated(topomesh%cell%vertP)) then 
            deallocate(topomesh%cell%vertP)
        end if 
        if (allocated(topomesh%cell%face)) then 
            deallocate(topomesh%cell%face)
        end if 
        if (allocated(topomesh%cell%faceP)) then 
            deallocate(topomesh%cell%faceP)
        end if 

        topomesh%cell%ntot = cc
        topomesh%cell%ID = [(i, i = 1, cc)]
        if(allocated(topomesh%cell%flags)) then 
            deallocate(topomesh%cell%flags)
        end if 
        allocate(topomesh%cell%flags(cc))
        topomesh%cell%flags = 0_I8
        if (allocated(topomesh%cell%tube)) then 
            deallocate(topomesh%cell%tube)
        end if
        allocate(topomesh%cell%tube(cc))
        topomesh%cell%tube = 0_I8
        allocate(topomesh%cell%vertP(cc, 2))
        allocate(topomesh%cell%faceP(cc, 2))
        topomesh%cell%vertP(1, 1) = 1
        topomesh%cell%faceP(1, 1) = 1
        do i = 1, cc
            topomesh%cell%vertP(i, 2) = cellvert(i)%size()
            topomesh%cell%faceP(i, 2) = cellface(i)%size()
        end do 
        allocate(topomesh%cell%vert(sum(topomesh%cell%vertP(:, 2))))
        allocate(topomesh%cell%face(sum(topomesh%cell%faceP(:, 2))))
        do i = 2, cc
            topomesh%cell%vertP(i, 1) = topomesh%cell%vertP(i-1, 1) + &
                topomesh%cell%vertP(i-1, 2)
            topomesh%cell%faceP(i, 1) = topomesh%cell%faceP(i-1, 1) + &
                topomesh%cell%faceP(i-1, 2)
        end do 
        do i = 1, cc
            topomesh%cell%vert(&
                topomesh%cell%vertP(i, 1):topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1) = & 
                cellvert(i)%Get()
            topomesh%cell%face(&
                topomesh%cell%faceP(i, 1):topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1) = & 
                cellface(i)%Get()
        end do 
        topomesh%cell%ntotv = size(topomesh%cell%vert)
        topomesh%cell%ntotf = size(topomesh%cell%face)

        ! Housekeeping
        !=============
        end associate
        
    end subroutine 

    ! Data addition
    subroutine AddTopologicalMeshData(topomesh)

        ! Description
        !============
        ! Add data of the vertices, faces, and cells (if they are present) based on
        ! their type. The following fields are updated:
        ! - vert:
        !   * BV, IV: logical (nv-by-1) for boundary or internal vertex. Boundary
        !   vertices are those vertices with type 4, 5 or 6, all the rest are
        !   internal. (i.e. tangency points and regular boundary points)
        ! - face:
        !   * BF, IF: logical (nf-by-1) for boundary or internal face. Boundary
        !   faces are those faces with type 3 or 7, all the rest are internal (unless
        !   cells have already been computed and we can check the amount of 
        !   cells per face)

        ! Note: it is assumed that all necessary contours/faces have 
        ! been added already
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Loop
        integer(I8)                             :: i 

        ! Vertices
        !=========
        ! Logicals
        topomesh%vert%BV = (topomesh%vert%type == TMvertextp1ID) .or. &
            (topomesh%vert%type == TMvertextp2ID) .or. &
            (topomesh%vert%type == TMvertexbndID)

        ! Faces
        !======
        ! Logicals
        topomesh%face%BF = (topomesh%face%type == TMfacebndID ) ! just for init
        topomesh%face%BF = .false. 
        do i = 1, size(TMfaceBFID)
            topomesh%face%BF = topomesh%face%BF .or. &
                topomesh%face%type == TMfaceBFID(i)
        end do 
        

        ! Check if cells are present
        if (allocated(topomesh%face%cellP)) then 
            if (size(topomesh%face%cellP, 1) == topomesh%face%ntot) then 
                topomesh%face%BF = topomesh%face%cellP(:, 2) == 1
            end if 
        end if 

        ! Labels
        if (allocated(topomesh%face%label)) then 
            deallocate(topomesh%face%label)
        end if 
        allocate(topomesh%face%label(topomesh%face%ntot))
        do i = 1, topomesh%face%ntot
            topomesh%face%label(i) = i
        end do 

    end subroutine

    ! Flux tube addition
    subroutine AddTopologicalMeshTubes(topomesh, types)

        ! Description
        !============
        ! This routine attempts to extract 'tubes' of cells (e.g. flux tubes for
        ! fusion applications). These tubes are sequences of cells that have
        ! exactly two faces of a certain topological type (or set of types). Tubes
        ! are either closed (starting and ending at the same cell) or open
        ! (starting and ending at boundary faces). If tubes can't be constructed,
        ! the output is the empty struct. If some cells remain (i.e. some tubes
        ! could not be constructed), a message is shown. 

        ! Note: 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: types(:)

        ! Auxiliary
        integer(I8)                             :: thistf, thisc, &
            nextc, tc, ne, si1, si2, cpID
        integer(I8), allocatable, dimension(:)  :: tf, thesecells, &
            tfc, thiscf, tv, tfbnd, sortindex, sortedbndfaces, tf1, tf2, &
            tv1, tv2, polygonID
        integer(I8), allocatable, dimension(:, :)   :: bndfacevert, tfv
        logical                                 :: isclosed, addtube, &
            singlefacetube
        logical, allocatable, dimension(:)      :: hasfacetype, &
            hasfaces, iscellnotfound, ispolygonstart, isbranchingpolygon
        type(IntegerDynamicArrayUDT)            :: tubef, tubec, &
            ntubef, ntubec, isclosedtube, alltubef, alltubec, &
            tubebndf, alltubebndf1, alltubebndf2, ntubebndf1, ntubebndf2, &
            alltubebndv1, alltubebndv2, ntubebndv1, ntubebndv2, &
            temptf1, temptf2

        ! Loop
        integer(I8)                             :: i, k, ct

        ! Initialize
        !===========
        ! Unpack
        associate(&
            face            => topomesh%face,   &
            cell            => topomesh%cell    &
            )


        ! Get faces & cells
        !==================
        ! Check faces
        allocate(hasfacetype(face%ntot))
        hasfacetype = .false. 
        do i = 1, size(types)
            hasfacetype = hasfacetype .or. (face%type == types(i));
        end do

        ! Check cells
        allocate(hasfaces(cell%ntot))
        hasfaces = .false. 
        do i = 1, cell%ntot
            ! Get face cells
            tf = GetTMCellFace(cell, i)
            if (count(hasfacetype(tf)) == 2) then 
                hasfaces(i) = .true.
            end if 
        end do 

        ! Display warning if some cells were not found
        if (any(.not. hasfaces)) then 
            allocate(thesecells(count(.not. hasfaces)))
            thesecells = pack([(k, k = 1, cell%ntot)], .not. hasfaces)
            print *, 'The following cells do not have exactly two ' // & 
                'tube faces and will not be part of a tube: ',  thesecells
            print *, 'Vertices of these cells:'
            do i = 1, size(thesecells)
                tv = GetTMCellVert(topomesh%cell, thesecells(i))
                print *, 'vertices: ', tv 
            end do 
            call WriteTopologicalMesh(topomesh, 'topomesh_tube', .false.)
            deallocate(thesecells)
        end if 

        ! Construct tubes
        !================
        ! Initialize
        iscellnotfound = hasfaces
        ct = 0
        ntubef          = ConstructIntegerDynamicArray()
        ntubec          = ConstructIntegerDynamicArray()
        alltubef        = ConstructIntegerDynamicArray()
        alltubebndf1    = ConstructIntegerDynamicArray()
        alltubebndf2    = ConstructIntegerDynamicArray()
        alltubebndv1    = ConstructIntegerDynamicArray()
        alltubebndv2    = ConstructIntegerDynamicArray()
        ntubebndv1      = ConstructIntegerDynamicArray()
        ntubebndv2      = ConstructIntegerDynamicArray()
        ntubebndf1      = ConstructIntegerDynamicArray()
        ntubebndf2      = ConstructIntegerDynamicArray()
        alltubec        = ConstructIntegerDynamicArray()
        isclosedtube    = ConstructIntegerDynamicArray()

        ! Loop
        do while (.true.)
            
            ! Take the first cell
            tc = findloc(iscellnotfound, .true., 1)
            if (tc == 0) then 
                ! all found, exit
                exit 
            end if 
            iscellnotfound(tc) = .false.
            
            ! Initialize a new tube
            tubef           = ConstructIntegerDynamicArray()
            tubec           = ConstructIntegerDynamicArray()
            tubebndf        = ConstructIntegerDynamicArray()
            
            isclosed        = .false. 
            singlefacetube  = .false.
            addtube         = .true.

            ! Append the cell
            call tubec%Append(tc)
            
            ! Get the faces
            tf = GetTMCellFace(cell, tc)

            ! Get the bounding faces
            allocate(tfbnd(count(.not. hasfacetype(tf))))
            tfbnd = pack(tf, .not. hasfacetype(tf))
            call tubebndf%Append(tfbnd)
            deallocate(tfbnd)
            
            ! Get only faces that have the type
            tf = pack(tf, hasfacetype(tf))
            
            ! Sanity check
            if (size(tf) /= 2) then 
                call gdErrorHandler('AddTopologicalMeshTubes: did not ' // & 
                    'find two faces for cell which should have two ' // & 
                    'faces - this is a bug')
            end if 

            ! Check if both faces are the same - then don't trace, single cell tube
            if (tf(1) == tf(2)) then 
                isclosed = .true.
                singlefacetube = .true.
            end if 
            
            ! Take the first face and trace
            thistf = tf(1)
            thisc = tc
            
            ! Add if it is a single cell tube
            if (singlefacetube) then 
                call tubef%Append(tf)
            end if 
            
            ! Loop until all flux tube cells found from one side
            do while (.true. .and. (.not. singlefacetube)) 
                ! Add the current face before already found faces (to keep correct
                ! order)
                call tubef%Insert(thistf, 1)
                
                ! Is the face a boundary face? 
                if (face%BF(thistf)) then 
                    ! Exit
                    exit 
                end if 
                
                ! Get face cells (should be exactly two)
                tfc = GetTMFaceCell(face, thistf)
                if (size(tfc) /= 2 ) then 
                    call gdErrorHandler('Did not find two cells for an ' // & 
                        'inner face - this is a bug')
                end if 
                
                ! Check which cell is next
                if (tfc(1) /= thisc) then 
                    nextc = tfc(1)
                elseif (tfc(2) /= thisc) then 
                    nextc = tfc(2) 
                else 
                    call gdErrorHandler('Did not find next cell - this ' // &
                        'is a bug')
                end if
                
                ! Check if this cell is marked, otherwise we need to stop here and
                ! not add the flux tube
                if (((.not. hasfaces(nextc)) .or. (.not. iscellnotfound(nextc))) &
                    .and. (nextc /= tc))   then 
                    print *, 'AddTopologicalMeshTubes: could not extract a tube'
                    addtube = .false.
                    exit 
                end if 
                
                ! If we got here, we can safely add the cell and find the next face
                thisc = nextc
                if (nextc /= tc) then 
                    call tubec%Insert(thisc, 1)
                end if 
                iscellnotfound(thisc) = .false.
                
                ! Get cell faces
                thiscf = GetTMCellFace(cell, thisc)

                ! Add bounding faces
                allocate(tfbnd(count(.not. hasfacetype(thiscf))))
                tfbnd = pack(thiscf, .not. hasfacetype(thiscf))
                call tubebndf%Append(tfbnd)
                deallocate(tfbnd)

                ! Find next face
                thiscf = pack(thiscf, hasfacetype(thiscf) .and. thiscf /= thistf)
                
                ! Sanity check
                if (size(thiscf) > 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: too ' // & 
                        'many faces found when trying to find next face')
                elseif (size(thiscf) < 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: could ' // & 
                        'not find next face')
                end if 
                thistf = thiscf(1)
                
                ! Is the face the same face as we started from? (closed tube)
                if (thistf == tf(1)) then 
                    isclosed = .true.
                    call tubef%Append(tf(2)) ! add at end

                    ! Exit
                    exit 
                end if 
            end do 

            ! Check if we should cycle
            if (.not. addtube) then 
                cycle 
            end if 
            
            ! Take the second face and trace, unless we can't add the tube or the
            ! tube is closed
            thistf = tf(2)
            thisc = tc
            do while ((.true.) .and. (.not. isclosed))  
                ! Add the current face after already found faces (to keep correct
                ! order)
                call tubef%Append(thistf)

                ! Is the face a boundary face? 
                if (face%BF(thistf)) then 
                    ! Exit
                    exit 
                end if 

                ! Get face cells (should be exactly two)
                tfc = GetTMFaceCell(face, thistf)
                if (size(tfc) /= 2 ) then 
                    call gdErrorHandler('Did not find two cells for an ' // & 
                        'inner face - this is a bug')
                end if 
                
                ! Check which cell is next
                if (tfc(1) /= thisc) then 
                    nextc = tfc(1)
                elseif (tfc(2) /= thisc) then 
                    nextc = tfc(2) 
                else 
                    call gdErrorHandler('Did not find next cell - this ' // &
                        'is a bug')
                end if
                
                ! Check if this cell is marked, otherwise we need to stop here and
                ! not add the flux tube
                if (((.not. hasfaces(nextc)) .or. (.not. iscellnotfound(nextc))) &
                    .and. (nextc /= tc))   then 
                    print *, 'AddTopologicalMeshTubes: could not extract a tube'
                    addtube = .false.
                    exit 
                end if 
                
                ! If we got here, we can safely add the cell and find the next face
                thisc = nextc
                if (nextc /= tc) then 
                    call tubec%Append(thisc)
                end if 
                iscellnotfound(thisc) = .false.
                
                ! Get cell faces
                thiscf = GetTMCellFace(cell, thisc)

                ! Add bounding faces
                allocate(tfbnd(count(.not. hasfacetype(thiscf))))
                tfbnd = pack(thiscf, .not. hasfacetype(thiscf))
                call tubebndf%Append(tfbnd)
                deallocate(tfbnd)

                ! Find the next face
                thiscf = pack(thiscf, hasfacetype(thiscf) .and. thiscf /= thistf)
                
                ! Sanity check
                if (size(thiscf) > 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: too ' // & 
                        'many faces found when trying to find next face')
                elseif (size(thiscf) < 1) then 
                    call gdErrorHandler('AddTopologicalMeshTubes: could ' // & 
                        'not find next face')
                end if 
                thistf = thiscf(1)
                
            end do 
            
            ! Add data to the tube if possible
            if (addtube) then 
                ! Update counter
                ct = ct + 1

                ! Add number of faces and cells
                call ntubef%Append(tubef%Size())
                call ntubec%Append(tubec%Size())
                if (isclosed) then 
                    call isclosedtube%Append(1_I8)
                else 
                    call isclosedtube%Append(0_I8)
                end if 

                ! Add faces to total faces etc
                call alltubef%Append(tubef%Get())
                call alltubec%Append(tubec%Get())

                ! Get 'polygons' from bounding faces
                ne = tubebndf%Size() 
                bndfacevert = face%vert(tubebndf%Get(), :)
                allocate(sortindex(ne), ispolygonstart(ne), isbranchingpolygon(ne))
                call SortPolygonEdges(bndfacevert, ne, sortindex, &
                    ispolygonstart, isbranchingpolygon, polygonID)
                bndfacevert(:, 1) = bndfacevert(sortindex, 1)
                bndfacevert(:, 2) = bndfacevert(sortindex, 2)
                sortedbndfaces = tubebndf%Get(sortindex)
                
                
                ! Checks
                if (any(isbranchingpolygon) .and. ne > 0) then 
                    ! Branching polygons, need to see which ones belong
                    ! together (max 2)
                    if (maxval(polygonID) > 2) then 
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube has more than two branching polygons as ' // & 
                            'boundary, unexpected')
                    end if

                    ! Extract polygon faces
                    cpID = 0
                    temptf1 = ConstructIntegerDynamicArray()
                    temptf2 = ConstructIntegerDynamicArray()
                    do i = 1, ne
                        ! Get the polygon ID if it is a start
                        if (ispolygonstart(i)) then 
                            cpID = polygonID(i)
                        end if 

                        ! Add face
                        if (cpID == 1_I8) then 
                            call temptf1%Append(sortedbndfaces(i))
                        elseif (cpID == 2_I8) then 
                            call temptf2%Append(sortedbndfaces(i))
                        else
                            ! We shouldn't get here
                            call gdErrorHandler('Unexpected error')
                        end if 
                    end do 

                    ! At least the first polygon should exist
                    call alltubebndf1%Append(temptf1%Get())
                    call ntubebndf1%Append(temptf1%Size())

                    ! Add first set
                    tf1 = temptf1%Get()

                    ! Extract vertices of first set - unsorted here...
                    allocate(tv1(size(tf1)+1))
                    bndfacevert = face%vert(tf1, :)
                    call Unique([bndfacevert(:, 1), bndfacevert(:, 2)], tv1)
                    call alltubebndv1%Append(tv1)
                    call ntubebndv1%Append(size(tv1))
                    deallocate(tv1)

                    
                    ! Check if the second set exists
                    if (temptf2%Size() > 0) then 
                        ! Extract faces
                        call alltubebndf2%Append(temptf2%Get())
                        call ntubebndf2%Append(temptf2%Size())

                        ! Add second set
                        tf2 = temptf2%Get()

                        ! Extract vertices of second set 
                        allocate(tv2(size(tf2)+1))
                        bndfacevert = face%vert(tf2, :)
                        call Unique([bndfacevert(:, 1), bndfacevert(:, 2)], tv2)
                        call alltubebndv2%Append(tv2)
                        call ntubebndv2%Append(size(tv2))
                        deallocate(tv2)
                    else
                        ! Append a zero
                        call ntubebndf2%Append(0_I8)

                        ! Extract vertices
                        call ntubebndv2%Append(1_I8)
                        tf = tubef%Get()
                        tfv = face%vert(tf, :) 
                        if (all(tfv(1, 1) == tfv(:, 1) .or. tfv(1, 1) == tfv(:, 2))) then 
                            call alltubebndv2%Append(tfv(1, 1))
                        elseif (all(tfv(1, 2) == tfv(:, 1) .or. tfv(1, 2) == tfv(:, 2))) then
                            call alltubebndv2%Append(tfv(1, 2))
                        else
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshTubes: ' // &
                                'found only one set of faces, but could not ' // & 
                                'find a common point of all tube faces as ' // & 
                                'other boundary side. Unexpected.')
                        end if 
                    end if 

                else 
                    ! Simple polygons, do sanity checks
                    if ((count(ispolygonstart) < 1) .and. (ne > 0)) then 
                        ! May be supported in the future
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube does not have at least one aligned face, ' // & 
                            'not supported.')
                    elseif ((count(ispolygonstart) > 2) .and. (ne > 0)) then 
                        ! Weird 
                        call gdErrorHandler('AddTopologicalMeshTubes: ' // & 
                            'tube has more than two aligned face boundaries, ' // & 
                            'unexpected.')
                    end if 

                    ! Extract boundary faces and vertices
                    if (count(ispolygonstart) == 2) then 
                        ! Get starting indices (normally first one is always true)
                        si1 = 1
                        si2 = findloc(ispolygonstart, .true., 1, back=.true.)

                        ! Add first set
                        tf1 = sortedbndfaces(si1:si2-1)
                        call alltubebndf1%Append(tf1)
                        call ntubebndf1%Append(size(tf1))

                        ! Add second set
                        tf2 = sortedbndfaces(si2:ne)
                        call alltubebndf2%Append(tf2)
                        call ntubebndf2%Append(size(tf2))

                        ! Extract vertices of first set 
                        allocate(tv1(size(tf1)+1))
                        call ExtractPolygonVertices(bndfacevert(1:si2-1, 1:2), &
                            size(tf1), tv1)
                        call alltubebndv1%Append(tv1)
                        call ntubebndv1%Append(size(tv1))
                        deallocate(tv1)

                        ! Extract vertices of second set 
                        allocate(tv2(size(tf2)+1))
                        call ExtractPolygonVertices(bndfacevert(si2:ne, 1:2), &
                            size(tf2), tv2)
                        call alltubebndv2%Append(tv2)
                        call ntubebndv2%Append(size(tv2))
                        deallocate(tv2)

                    elseif (count(ispolygonstart) == 1) then 
                        ! Only one segment, need to search for tangency point
                        ! vertex (should be one vertex in common with all 
                        ! tube radial faces)

                        ! Add first set/
                        tf1 = sortedbndfaces
                        call alltubebndf1%Append(tf1)
                        call ntubebndf1%Append(size(tf1))

                        ! 'Add' second set
                        call ntubebndf2%Append(0_I8)

                        ! Extract vertices of first set 
                        allocate(tv1(size(tf1)+1))
                        call ExtractPolygonVertices(bndfacevert, &
                            size(tf1), tv1)
                        call alltubebndv1%Append(tv1)
                        call ntubebndv1%Append(size(tv1))
                        
                        ! Extract vertices of second set 
                        call ntubebndv2%Append(1_I8)
                        tf = tubef%Get()
                        tfv = face%vert(tf, :) 
                        if (any(tfv(1, 1) == tv1)) then 
                            call alltubebndv2%Append(tfv(1, 2))
                        elseif (any(tfv(1, 2) == tv1)) then
                            call alltubebndv2%Append(tfv(1, 1))
                        else
                            ! This shouldn't happen
                            call gdErrorHandler('AddTopologicalMeshTubes: ' // &
                                'found only one set of faces, but could not ' // & 
                                'find a common point of all tube faces as ' // & 
                                'other boundary side. Unexpected.')
                        end if 
                        deallocate(tv1)
                    end if 
                end if 

                ! Housekeeping
                deallocate(sortindex, ispolygonstart, isbranchingpolygon)
            end if 
            
        end do 

        ! Construct
        !==========
        ! Ensure deallocation
        call topomesh%tube%Deallocate()

        ! Reallocate
        call topomesh%tube%Initialize(ct, alltubef%Size(), alltubec%Size(), &
            alltubebndf1%Size(), alltubebndf2%Size(), alltubebndv1%Size(), &
            alltubebndv2%Size())

        ! Set fields
        topomesh%tube%face      = alltubef%Get()
        topomesh%tube%cell      = alltubec%Get()
        topomesh%tube%bndf1     = alltubebndf1%Get()
        topomesh%tube%bndf2     = alltubebndf2%Get()
        topomesh%tube%bndv1     = alltubebndv1%Get()
        topomesh%tube%bndv2     = alltubebndv2%Get()
        topomesh%tube%isclosed  = isclosedtube%Get() == 1_I8

        topomesh%tube%faceP(:, 2)   = ntubef%Get()
        topomesh%tube%cellp(:, 2)   = ntubec%Get()
        topomesh%tube%bndf1P(:, 2)  = ntubebndf1%Get()
        topomesh%tube%bndf2P(:, 2)  = ntubebndf2%Get()
        topomesh%tube%bndv1P(:, 2)  = ntubebndv1%Get()
        topomesh%tube%bndv2P(:, 2)  = ntubebndv2%Get()
        if (topomesh%tube%ntot > 0) then 
            topomesh%tube%faceP(1, 1)   = 1
            topomesh%tube%cellP(1, 1)   = 1
            topomesh%tube%bndf1P(1, 1)  = 1
            topomesh%tube%bndf2P(1, 1)  = 1
            topomesh%tube%bndv1P(1, 1)  = 1
            topomesh%tube%bndv2P(1, 1)  = 1
        end if 

        do i = 2, topomesh%tube%ntot 
            topomesh%tube%faceP(i, 1) = topomesh%tube%faceP(i-1, 1) + &     
                topomesh%tube%faceP(i-1, 2)
            topomesh%tube%cellP(i, 1) = topomesh%tube%cellP(i-1, 1) + &     
                topomesh%tube%cellP(i-1, 2)    
            topomesh%tube%bndf1P(i, 1) = topomesh%tube%bndf1P(i-1, 1) + &     
                topomesh%tube%bndf1P(i-1, 2)
            topomesh%tube%bndf2P(i, 1) = topomesh%tube%bndf2P(i-1, 1) + &     
                topomesh%tube%bndf2P(i-1, 2)    
            topomesh%tube%bndv1P(i, 1) = topomesh%tube%bndv1P(i-1, 1) + &     
                topomesh%tube%bndv1P(i-1, 2)
            topomesh%tube%bndv2P(i, 1) = topomesh%tube%bndv2P(i-1, 1) + &     
                topomesh%tube%bndv2P(i-1, 2)    
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Flux tube data addition
    subroutine AddTopologicalMeshTubeData(topomesh)

        ! Description
        !============
        ! This routine adds additional data to the topological mesh 
        ! tubes, such as which cells are in a tube (cell%tube) and which
        ! neighbours a tube has. Also, a directional graph
        ! representation is built, where sequence is from high psi 
        ! values to low psi values. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(inout)       :: topomesh 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: ttc, tf1, tf2, tc, &
            tf1u, tf2u, ev1, ev2, v, tnb1, tnb2, tv1, tv2
        real(R8), allocatable, dimension(:)     :: tpsi1, tpsi2

        type(IntegerDynamicArrayUDT)            :: ftneig1, ftneig2

        ! Loop
        integer(I8)                             :: i, j, k, nnb1, nnb2

        ! Initialize
        !===========
        ! Associate for ease
        associate(&
            tube        => topomesh%tube, &
            cell        => topomesh%cell, &
            face        => topomesh%face, &
            vert        => topomesh%vert)

        ! Allocate 
        if (allocated(tube%ftneig1)) then 
            deallocate(tube%ftneig1)
        end if 
        if (allocated(tube%ftneig2)) then 
            deallocate(tube%ftneig2)
        end if 
        if (allocated(tube%ftneig1P)) then 
            deallocate(tube%ftneig1P)
        end if 
        if (allocated(tube%ftneig2P)) then 
            deallocate(tube%ftneig2P)
        end if 
        if (allocated(cell%tube)) then 
            deallocate(cell%tube)
        end if 
        allocate(cell%tube(cell%ntot))
        allocate(tube%ftneig1P(tube%ntot, 2), tube%ftneig2P(tube%ntot, 2))

        ! Initialize
        ftneig1     = ConstructIntegerDynamicArray()
        ftneig2     = ConstructIntegerDynamicArray()
        cell%tube       = 0
        tube%ftneig1P = 0
        tube%ftneig2P = 0
        if (tube%ntot > 1) then 
            tube%ftneig1P(1, 1) = 1_I8
            tube%ftneig2P(1, 1) = 1_I8
        end if  
        
        ! Cell tubes
        !===========
        do i = 1, tube%ntot
            ! Get tube cells
            ttc = tube%GetCell(i)

            ! Check
            if (any(cell%tube(ttc) /= 0)) then 
                print *, 'AddTopologicalMeshTubeData: some cells belong ' // &
                    'to multiple tubes, unexpected. Overwriting and continuing...'
            end if 

            ! Set cell tube
            cell%tube(ttc) = i
        end do 

        ! Tube neighbours
        !================
        ! Initialize
        nnb1 = 0
        nnb2 = 0

        ! Loop
        do i = 1, tube%ntot
            ! Get tube faces from first side
            tf1 = tube%GetBndFace(i, 1_I8)
            if (size(tf1) > 0) then 
                ! Loop over all face
                do j = 1, size(tf1)
                    ! Get the face cells
                    tc = face%GetCell(tf1(j))

                    ! Check if any cell tube IDs are not equal to the 
                    ! current ID or zero. If so, add neighbour
                    do k = 1, size(tc)
                        if (all(cell%tube(tc(k)) /= [i, 0])) then 
                            call ftneig1%Append(cell%tube(tc(k)))
                            tube%ftneig1P(i, 2) = tube%ftneig1P(i, 2) + 1
                        end if
                    end do
                end do
            end if 

            ! Get tube faces from second side
            tf2 = tube%GetBndFace(i, 2_I8)
            if (size(tf2) > 0) then 
                ! Loop over all face
                do j = 1, size(tf2)
                    ! Get the face cells
                    tc = face%GetCell(tf2(j))

                    ! Check if any cell tube IDs are not equal to the 
                    ! current ID. If so, add neighbour
                    do k = 1, size(tc)
                        if (all(cell%tube(tc(k)) /= [i, 0])) then 
                            call ftneig2%Append(cell%tube(tc(k)))
                            tube%ftneig2P(i, 2) = tube%ftneig2P(i, 2) + 1
                        end if
                    end do
                end do
            end if 
        end do 

        ! Build (preliminary)
        tube%ftneig1 = ftneig1%Get()
        tube%ftneig2 = ftneig2%Get()
        do i = 2, tube%ntot
            tube%ftneig1P(i, 1) = tube%ftneig1P(i-1, 1) + tube%ftneig1P(i-1, 2)
            tube%ftneig2P(i, 1) = tube%ftneig2P(i-1, 1) + tube%ftneig2P(i-1, 2)
        end do 

        ! Reduce to unique set for each tube
        ftneig1 = ConstructIntegerDynamicArray()
        ftneig2 = ConstructIntegerDynamicArray()
        do i = 1, tube%ntot
            ! Get tube neighbours (this is possible now)
            tf1 = tube%GetNeig(i, 1_I8)
            tf2 = tube%GetNeig(i, 2_I8)

            ! Determine unique set
            call Unique(tf1, tf1u)
            call Unique(tf2, tf2u)

            ! Append
            call ftneig1%Append(tf1u)
            call ftneig2%Append(tf2u)

            ! Set pointer counter, but only of this tube
            tube%ftneig1P(i, 2) = size(tf1u)
            tube%ftneig2P(i, 2) = size(tf2u)
        end do 

        ! Build (final)
        tube%ftneig1 = ftneig1%Get()
        tube%ftneig2 = ftneig2%Get()
        do i = 2, tube%ntot
            tube%ftneig1P(i, 1) = tube%ftneig1P(i-1, 1) + tube%ftneig1P(i-1, 2)
            tube%ftneig2P(i, 1) = tube%ftneig2P(i-1, 1) + tube%ftneig2P(i-1, 2)
        end do 

        ! Graph
        !======
        ! Vertices 
        v = [(k, k = 1, tube%ntot)] ! simple 

        ! Edges
        allocate(ev1(0), ev2(0))
        do i = 1, tube%ntot 
            ! Get all neighbours of this tube
            tnb1 = GetTMTubeNeig(tube, i, 1_I8)
            tnb2 = GetTMTubeNeig(tube, i, 2_I8)

            ! Get all faces
            tf1 = GetTMTubeBndFace(tube, i, 1_I8)
            tf2 = GetTMTubeBndFace(tube, i, 2_I8)

            ! Get all vertices
            tv1 = GetTMTubeBndVert(tube, i, 1_I8)
            tv2 = GetTMTubeBndVert(tube, i, 2_I8)
            
            if (any(face%fsID(tf1) == 0) .or. any(face%fsID(tf2) == 0)) then 
                print *, 'something weird'
                call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
            end if 

            ! Get all flux values (include vertex values for tangency point cases)
            tpsi1 = [topomesh%fsfval%Get(face%fsID(tf1)), topomesh%vert%fval(tv1)]
            tpsi2 = [topomesh%fsfval%Get(face%fsID(tf2)), topomesh%vert%fval(tv2)]

            ! Sanity checks
            if (all(minval(tpsi1) > tpsi2)) then 
                ! First boundary is high field boundary -> first neighbours
                ! should be high field neighbours, so edges go from these 
                ! neighbours to the current tube. Second neighbours are
                ! low field, so edges go from this tube to nb2 tubes
                ev1 = [ev1, tnb1]
                ev2 = [ev2, spread(i, 1, size(tnb1))]
                
                ev1 = [ev1, spread(i, 1, size(tnb2))]
                ev2 = [ev2, tnb2]

                ! Also set hfside
                tube%hfside(i) = 1

            elseif (all(minval(tpsi2) > tpsi1)) then 
                ! Second boundary is high field boundary 
                ev1 = [ev1, tnb2]
                ev2 = [ev2, spread(i, 1, size(tnb2))]
                
                ev1 = [ev1, spread(i, 1, size(tnb1))]
                ev2 = [ev2, tnb1]

                ! Also set hfside
                tube%hfside(i) = 2
            else
                ! May originate due to merging. In that case, compare
                ! extremal values
                call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
                print *, 'vertices: ', [tv1, tv2]
                print *, 'AddTopologicalMeshTubeData: ' // & 
                    'boundaries of tube have different psi values and ' // & 
                    'overlap - trying to distinguish based on maximal and ' // & 
                    'minimal value. This may happen during tube merging...'
                if (all(maxval(tpsi1) > tpsi2)) then 
                    ! First boundary is high field boundary 
                    ev1 = [ev1, tnb1]
                    ev2 = [ev2, spread(i, 1, size(tnb1))]
                    
                    ev1 = [ev1, spread(i, 1, size(tnb2))]
                    ev2 = [ev2, tnb2]

                    ! Also set hfside
                    tube%hfside(i) = 1
                else
                    ! Second boundary is high field boundary 
                    ev1 = [ev1, tnb2]
                    ev2 = [ev2, spread(i, 1, size(tnb2))]
                    
                    ev1 = [ev1, spread(i, 1, size(tnb1))]
                    ev2 = [ev2, tnb1]

                    ! Also set hfside
                    tube%hfside(i) = 2
                end if 

                

            end if 

        end do 

        ! Construct
        call tube%graph%Construct(ev1, ev2, v)

        ! Condense to remove duplicated edges
        call tube%graph%Condense()

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Interconnection
    subroutine AddTopologicalMeshInterconnectionData(topomesh)

        ! Description
        !============
        ! Add additional interconnection data to the topological mesh. 
        ! Similar to AddGridInterconnections. Some fields may be 
        ! constructed again that already existed (e.g. cell faces etc), 
        ! but we keep it in here to have a general interconnections 
        ! routine.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 

        ! Auxiliary
        integer(I8), allocatable                :: tempvcells(:, :)
        integer(I8), allocatable, dimension(:)  :: ncpv, vcount, tv

        ! Loop
        integer(I8)                             :: i, j

        ! Add face cells
        !===============
        ! Because separate routine
        call AddTopologicalMeshFaceCells(topomesh)
         
        ! Unpack & initialize
        !====================
        ! Data structures
        associate(&
            nv      => topomesh%vert%ntot,  &
            nf      => topomesh%face%ntot,  &
            nc      => topomesh%cell%ntot,  &
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell    &
            )
        
        ! Checks
        if (any((f%vert(:, 1) == 0) .or. (f%vert(:, 2) == 0))) then 
            call gdErrorHandler('AddTopologicalMeshInterconnectionData : ' // & 
                'Some vertex indices are zero in faces, not supported')
        end if 
        if (allocated(v%cell)) then 
            deallocate(v%cell)
        end if 
        if (allocated(v%cellP)) then 
            deallocate(v%cellP)
        end if 
        
        ! Basic interconnections
        !=======================
        ! Initialize
        allocate(ncpv(v%ntot))
        ncpv = 0
        do i = 1, size(c%vert)
            ncpv(c%vert(i)) = ncpv(c%vert(i))+1
        end do
        allocate(v%cell(sum(ncpv)), v%cellP(v%ntot, 2))
        v%cellP(1, 1) = 1
        v%cellP(:, 2) = ncpv
        do i = 2, v%ntot
            v%cellP(i, 1) = v%cellP(i-1, 1) + v%cellP(i-1, 2)
        end do 
        
        ! Sanity check
        if (any(v%cellP(:, 2) < 1)) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_error', .false.)
            print *, findloc(v%cellP(:, 2), 0, 1)
            call gdErrorHandler('AddTopologicalMeshInterconnectionData: ' // & 
                'vertex without any cells detected. Check grid interconnectivity')
        end if 
        
        ! Construct vertex cells (cells of faces and faces of cells 
        ! already found in constructor phase)
        ! Note: we construct first
        ! temporary arrays (nv-by-ncpv, nf-by-2) that are afterwards converted to
        ! cell and cellP arrays. 
        allocate(tempvcells(v%ntot, maxval(ncpv)), vcount(v%ntot))
        tempvcells = 0
        vcount = 0
        do i = 1, nc
            ! Get vertices of cell
            tv = GetTMCellVert(c, i) ! there may be doubles in here!
            
            ! Add vertex cells
            do j = 1, size(tv)
                ! Update counter
                vcount(tv(j)) = vcount(tv(j))+1

                ! Add cell
                tempvcells(tv(j), vcount(tv(j))) = i 
            end do 
        end do 
        
        ! Construct 
        do i = 1, v%ntot
            v%cell(v%cellP(i, 1):v%cellP(i, 1)+v%cellP(i, 2)-1) = & 
                tempvcells(i, 1:vcount(i))
        end do 
        
        ! Add logicals
        f%BF = f%cellP(:, 2) == 1

        ! Reconstruct boundary vertices based on boundary faces
        if (allocated(v%BV)) then 
            deallocate(v%BV)
        end if
        allocate(v%BV(v%ntot))
        v%BV = .false. 
        do i = 1, f%ntot
            if (f%BF(i)) then 
                v%BV(f%vert(i, :)) = .true. 
            end if 
        end do 

        ! Reconstruct fsID of vertices to be compliant to faces (note: 
        ! type 1 tangency points get their own flux surface ID)

                    

        ! Tubes
        !======
        ! Flux tubes 
        call AddTopologicalMeshTubes(topomesh, [TMfaceradID, TMfacebndID])

        ! Additional tube interconnection data
        call AddTopologicalMeshTubeData(topomesh)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Vertex removal
    subroutine RemoveTopologicalMeshVertexLogical(topomesh, rmvert)

        ! Description
        !============
        ! Remove the vertices with IDs specified in 'IDs' from the topological
        ! mesh. Only vertex and face data is updated. Faces are not immediately
        ! removed, but the vertices of each face are updated (set to zero if
        ! vertex deleted, otherwise updated to new vertices). IDs is a 
        ! logical array of size 1:topomesh%vert%ntot

        ! Notes
        !======
        ! Note 1: Faces do not necessarily require start and end point yet (they
        ! can have zeros in face.vert)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, intent(in)                     :: rmvert(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: diffID, tv
        logical, allocatable, dimension(:)      :: keepvert

        ! Loop 
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Check
        if (size(rmvert) /= topomesh%vert%ntot) then 
            call gdErrorHandler('RemoveTopologicalMeshVertexLogical: ' // & 
                'rmvert has illegal size')
        end if 

        ! Initialize
        keepvert = .not. rmvert 

        ! Associate
        associate(&
            nv      => topomesh%vert%ntot,  &
            nf      => topomesh%face%ntot,  &
            nc      => topomesh%cell%ntot)

        ! Delete vertices
        !================
        ! Compute vertex ID adjuster
        allocate(diffID(nv))
        diffID = 0
        do i = 1, nv-1 
            if (rmvert(i)) then 
                diffID(i+1) = diffID(i) + 1
            else
                diffID(i+1) = diffID(i)
            end if
        end do

        ! Remove vertices
        topomesh%vert%ntot  = count(keepvert)
        topomesh%vert%ID    = topomesh%vert%ID - diffID 
        topomesh%vert%ID    = pack(topomesh%vert%ID, keepvert)
        topomesh%vert%x     = pack(topomesh%vert%x, keepvert)
        topomesh%vert%y     = pack(topomesh%vert%y, keepvert)
        topomesh%vert%fval  = pack(topomesh%vert%fval, keepvert)
        topomesh%vert%type  = pack(topomesh%vert%type, keepvert) 
        topomesh%vert%fsID  = pack(topomesh%vert%fsID, keepvert)

        ! Adjust faces
        !=============
        ! Loop over all faces
        do i = 1, nf
            ! First vertex
            if (topomesh%face%vert(i, 1) /= 0) then 
                if (rmvert(topomesh%face%vert(i, 1))) then 
                    topomesh%face%vert(i, 1) = 0
                else
                    topomesh%face%vert(i, 1) = topomesh%face%vert(i, 1) - diffID(topomesh%face%vert(i, 1))
                end if
            end if 
            if (topomesh%face%vert(i, 2) /= 0) then 
                if (rmvert(topomesh%face%vert(i, 2))) then 
                    topomesh%face%vert(i, 2) = 0
                else
                    topomesh%face%vert(i, 2) = topomesh%face%vert(i, 2) - diffID(topomesh%face%vert(i, 2))
                end if 
            end if
        end do

        ! Adjust cells
        !=============
        ! Loop over all cells
        do i = 1, nc 
            ! Get vertices
            tv = GetTMCellVert(topomesh%cell, i)
            
            ! Check
            do j = 1, size(tv)
                if (tv(j) /= 0) then
                    if (rmvert(tv(j))) then 
                        tv(j) = 0
                    else
                        tv(j) = tv(j) - diffID(tv(j))
                    end if 
                end if
            end do 
            !tv(delind(tv)) = 0;
            !tv(tv ~= 0) = tv(tv ~= 0) - diffID(tv(tv ~= 0));
            
            ! Reset
            topomesh%cell%vert(topomesh%cell%vertP(i, 1):topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1) = tv
        end do

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Face removal
    subroutine RemoveTopologicalMeshFaceLogical(topomesh, rmface)

        ! Description
        !============
        ! Remove the faces with IDs specified in rmface from the mesh. 
        ! IMPORTANT: no cell or vertex information is updated here! This should be
        ! recomputed afterwards. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        logical, intent(in)                 :: rmface(:)

        ! Auxiliary
        integer(I8), allocatable            :: tmp(:, :), diffIDf(:), &
            tcf(:)
        logical, allocatable, dimension(:)  :: keepface

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            nf  => topomesh%face%ntot,    &
            nc  => topomesh%cell%ntot)
        keepface = .not. rmface 

        ! Check
        if (size(rmface) /= nf) then 
            call gdErrorHandler('RemoveTopologicalMeshFaceLogical: ' // & 
                'illegal size of rmface')
        end if 
        if (any((topomesh%face%type == TMfacealbndID) .and. rmface)) then 
            call WriteTopologicalMesh(topomesh, 'topomesh_temp', .false.)
            print *, 'RemoveTopologicalMeshFaceLogical: removing aligned boundary faces'
        end if 

        ! Determine ID shift
        allocate(diffIDf(nf))
        diffIDf = 0
        do i = 1, nf-1 
            if (rmface(i)) then 
                diffIDf(i+1) = diffIDf(i) + 1
            else
                diffIDf(i+1) = diffIDf(i)
            end if
        end do

        ! Delete
        !=======
        ! Update IDs
        topomesh%face%ID = topomesh%face%ID - diffIDf

        ! Remove data
        topomesh%face%ntot = count(keepface)
        tmp = topomesh%face%vert 
        deallocate(topomesh%face%vert)
        allocate(topomesh%face%vert(count(keepface), 2))
        topomesh%face%vert(:, 1) = pack(tmp(:, 1), keepface)
        topomesh%face%vert(:, 2) = pack(tmp(:, 2), keepface)
        topomesh%face%x         = pack(topomesh%face%x, keepface)
        topomesh%face%y         = pack(topomesh%face%y, keepface)
        topomesh%face%fsID      = pack(topomesh%face%fsID, keepface)
        topomesh%face%type      = pack(topomesh%face%type, keepface)
        topomesh%face%ID        = pack(topomesh%face%ID, keepface)
        topomesh%face%pol       = pack(topomesh%face%pol, keepface)  
        
        if (allocated(topomesh%face%cell)) then 
            deallocate(topomesh%face%cell)
        end if 
        if (allocated(topomesh%face%cellP)) then 
            deallocate(topomesh%face%cellP)
        end if 

        ! Update cell face IDs
        do i = 1, nc
            ! Get faces
            tcf = GetTMCellFace(topomesh%cell, i)
            
            ! Reset faces (cells are not removed here!)
            topomesh%cell%face(topomesh%cell%faceP(i, 1):topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1) = tcf - diffIDf(tcf)
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Cell removal 
    subroutine RemoveTopologicalMeshCellLogical(topomesh, rmcell)

        ! Description
        !============
        ! Remove the cells with IDs specified in IDs from the mesh. 
        ! IMPORTANT: no face or vertex information is updated here! This should be
        ! recomputed afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        logical, intent(in)                     :: rmcell(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: indv, indf, vp, &
            fp
        logical, allocatable, dimension(:)      :: delv, delf 

        ! Loop
        integer(I8)                             :: i, k 

        ! Delete
        !=======
        ! Initialize
        allocate(delv(size(topomesh%cell%vert)), delf(size(topomesh%cell%face)))
        delv = .false. 
        delf = .false. 

        ! Mark vertices and faces for removal
        do i = 1, topomesh%cell%ntot 
            if (rmcell(i)) then 
                indv = [(k, k = topomesh%cell%vertP(i, 1), &
                    topomesh%cell%vertP(i, 1)+topomesh%cell%vertP(i, 2)-1)]
                indf = [(k, k = topomesh%cell%faceP(i, 1), &
                    topomesh%cell%faceP(i, 1)+topomesh%cell%faceP(i, 2)-1)]
                delf(indf) = .true.
                delv(indv) = .true.
            end if 
        end do 

        ! Remove
        topomesh%cell%ntot = count(.not. rmcell)
        topomesh%cell%vert = pack(topomesh%cell%vert, .not. delv)
        topomesh%cell%face = pack(topomesh%cell%face, .not. delf)
        vp = pack(topomesh%cell%vertP(:, 2), .not. rmcell)
        fp = pack(topomesh%cell%faceP(:, 2), .not. rmcell)
        deallocate(topomesh%cell%vertP, topomesh%cell%faceP)
        allocate(topomesh%cell%vertP(topomesh%cell%ntot, 2), &
            topomesh%cell%faceP(topomesh%cell%ntot, 2))
        topomesh%cell%vertP(:, 2) = vp 
        topomesh%cell%faceP(:, 2) = fp 
        topomesh%cell%vertP(1, 1) = 1
        topomesh%cell%faceP(1, 1) = 1
        do i = 2, topomesh%cell%ntot 
            topomesh%cell%vertP(i, 1) = topomesh%cell%vertP(i-1, 1) + &
                topomesh%cell%vertP(i-1, 2)
            topomesh%cell%faceP(i, 1) = topomesh%cell%faceP(i-1, 1) + &
                topomesh%cell%faceP(i-1, 2)
        end do
        topomesh%cell%ntotv = size(topomesh%cell%vert)
        topomesh%cell%ntotf = size(topomesh%cell%face)

    end subroutine

    ! Topological face extraction
    subroutine ExtractTopologicalFacesFromPolygon(pol, eID, sID, &
            xint, yint, v1, v2, xf, yf)

            ! Description
            !============
            ! This function extracts faces from polygons given the segment intersection
            ! and vertex IDs. If the points don't exactly coincide with one of the
            ! segment's vertices, the vertex is added. xint and yint should contain the
            ! coordinates of all intersections, of which the right ones can be queried
            ! with xint(IDs). It is assumed that all intersections are unique. 
            ! Furthermore, it is assumed that the segment intersection 
            ! indices, sID, are already properly sorted (this is not
            ! explicitly checked!)

            ! Notes
            !======
            ! Note 1: it is assumed that the intersections are unique.

            ! Note 2: close polygons are accounted for. It is assumed that the last
            ! point in the polygon.vert structure is the same as the first one. This is
            ! not added twice in facedata.x, y.

            ! Declare variables
            !==================
            ! Arguments
            type(PolygonUDT)                        :: pol 
            integer(I8), allocatable, intent(out)   :: v1(:), v2(:)
            integer(I8), intent(inout)              :: eID(:), sID(:)
            type(RealDynamicArrayUDT), allocatable, intent(out) :: xf(:), yf(:)
            real(R8), intent(in)                    :: xint(:), yint(:)

            ! Auxiliary
            integer(I8)                             :: nf, fc, &
                seID, eeID, ssID, esID 
            integer(I8), allocatable, dimension(:)  :: sortedsID(:), &
                sortedeID(:)
            real(R8)                                :: sx, sy, ex, ey
            real(R8), allocatable, dimension(:)     :: tx, ty

            ! Loop
            integer(I8)                             :: i, si, ei  

            ! Initialize
            !===========
            ! Check maximal number of faces
            nf = size(eID) + 1

            ! Allocate & initialize
            allocate(v1(nf), v2(nf), xf(nf), yf(nf))
            v1 = 0
            v2 = 0
            
            ! Face counter
            fc = 0 

            ! Check
            if (nf == 1) then 
                ! No intersections, just return polygon as is
                xf(1) = ConstructRealDynamicArray(pol%x(pol%vert))
                yf(1) = ConstructRealDynamicArray(pol%y(pol%vert))
                return 
            end if 

            ! Initialize
            do i = 1, nf 
                xf(i) = ConstructRealDynamicArray()
                yf(i) = ConstructRealDynamicArray()
            end do 

            ! Associate
            associate( &
                nsID    =>      size(sID),          & ! size(sID) == size(eID) normally
                neID    =>      size(eID),          &
                px      =>      pol%x(pol%vert),    &
                py      =>      pol%y(pol%vert)     &
                )

            ! Construct segments
            !===================
            sortedsID = sID ! assumed already sorted!
            sortedeID = eID

            ! Treat start and end segments
            if (pol%isclosed) then 
                ! Closed polygon
                if (( (sortedsID(1) == 1 .and. px(1) == xint(sortedeID(1)))) .and. & 
                    ( (sortedsID(nsID) == pol%ne .and. px(pol%ne+1) == xint(sortedeID(nsID))))) then 
                    ! Start and end exactly in start point, do nothing
                elseif (.not. ( (sortedsID(1) == 1 .and. px(1) == xint(sortedeID(1)))) .and. & 
                    .not. ( (sortedsID(nsID) == pol%ne .and. px(pol%ne+1) == xint(sortedeID(nsID))))) then 
                    ! Start and end segment without start and end vertex,
                    ! so merge
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = sortedeID(neID)
                    v2(fc) = sortedeID(1)

                    ! Add face coordinates (hedge for exact duplicates)
                    tx = [xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne)]
                    ty = [yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne)]
                    tx = [tx, px(1:sortedsID(1)), xint(sortedeID(1))]
                    ty = [ty, py(1:sortedsID(1)), yint(sortedeID(1))]
                    if ((tx(1) == tx(2)) .and. (ty(1) == ty(2))) then 
                        tx = tx(2:size(tx))
                        ty = ty(2:size(ty))
                    end if 
                    if ((tx(size(tx)-1) == tx(size(tx))) .and. &
                        (ty(size(ty)-1) == ty(size(ty)))) then 
                        tx = tx(1:size(tx)-1)
                        ty = ty(1:size(ty)-1)
                    end if 

                    call xf(fc)%Append(tx)
                    call yf(fc)%Append(ty)
                else
                    ! Something weird here: intersection in exactly first node, but
                    ! only at first and not at last segment. Likely because
                    ! intersection with last segment was not added. Throw error. 
                    call gdErrorHandler(&
                        'ExtractTopologicalFacesFromPolygon: ' // & 
                        'closed polygon intersects in start point ' // & 
                        'but not in end point - check input '// & 
                        '(end point intersection should be added separately)')
                end if 
            else 
                ! Open polygon
                if ((sortedsID(1) == 1) .and. (px(1) == xint(sortedeID(1))) & 
                    .and. (py(1) == yint(sortedeID(1)))) then 
                    ! Start veretx at first intersection, no starting segment
                    ! without starting vertex 
                else
                    ! Start at second intersection as indicated, add starting
                    ! segment 
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = 0 ! no start vertex
                    v2(fc) = sortedeID(1) ! end vertex 

                    ! Add face coordinates
                    tx = [px(1:sortedsID(1)), xint(sortedeID(1))]
                    ty = [py(1:sortedsID(1)), yint(sortedeID(1))]
                    if ((tx(size(tx)-1) == tx(size(tx))) .and. &
                        (ty(size(ty)-1) == ty(size(ty)))) then 
                        tx = tx(1:size(tx)-1)
                        ty = ty(1:size(ty)-1)
                    end if 
                    call xf(fc)%Append(tx) 
                    call yf(fc)%Append(ty)
                
                end if 
                if ((sortedsID(nsID) == pol%ne) .and. (px(pol%ne+1) == xint(sortedeID(nsID))) & 
                    .and. (py(pol%ne+1) == yint(sortedeID(nsID)))) then 
                    ! Start veretx at first intersection, no starting segment
                    ! without starting vertex 
                else
                    ! Start at second intersection as indicated, add starting
                    ! segment 
                    fc = fc + 1

                    ! Add vertices
                    v1(fc) = sortedeID(neID) ! start vertex
                    v2(fc) = 0 ! no end vertex 

                    ! Add face coordinates
                    tx = [xint(sortedeID(neID)), px(sortedsID(nsID)+1:pol%ne+1)]
                    ty = [yint(sortedeID(neID)), py(sortedsID(nsID)+1:pol%ne+1)]
                    if ((tx(1) == tx(2)) .and. (ty(1) == ty(2))) then 
                        tx = tx(2:size(tx))
                        ty = ty(2:size(ty))
                    end if 
                    call xf(fc)%Append(tx) 
                    call yf(fc)%Append(ty)
                
                end if 
            end if 

            ! Set start and end indices for other segments
            si = 1
            ei = neID-1

            ! Loop over remaining segments
            do i = si, ei 
                ! All segments should be nicely closed so no special 
                ! checks needed

                ! Get start and end points & segments
                seID = sortedeID(i)
                eeID = sortedeID(i+1)
                ssID = sortedsID(i)
                esID = sortedsID(i+1)
                sx = xint(seID)
                sy = yint(seID)
                ex = xint(eeID)
                ey = yint(eeID)

                ! Update counter
                fc = fc + 1

                ! Add vertices
                v1(fc) = seID 
                v2(fc) = eeID 

                ! Add coordinates
                tx = px(ssID+1:esID)
                ty = py(ssID+1:esID)

                if (.not. allocated(tx)) then 
                    allocate(tx(0), ty(0))
                end if 

                ! Check
                if (size(tx) > 0) then 
                    ! Check if we need to add start and end points
                    if ((tx(1) /= sx) .or. (ty(1) /= sy)) then 
                        tx = [sx, tx]
                        ty = [sy, ty]
                    end if 
                    if ((tx(size(tx)) /= ex) .or. (ty(size(ty)) /= ey)) then 
                        tx = [tx, ex]
                        ty = [ty, ey]
                    end if 
                else 
                    ! Only add start and end points
                    tx = [sx, ex]
                    ty = [sy, ey]
                end if 

                ! Add
                call xf(fc)%Append(tx)
                call yf(fc)%Append(ty)

            end do

            ! Check if we need to shrink (may happen for closed polygons)
            if (fc < nf) then 
                xf = xf(1:fc)
                yf = yf(1:fc)
                V1 = V1(1:fc)
                V2 = V2(1:fc)
            end if 

            ! Housekeeping
            !=============
            end associate





    end subroutine 

    ! Getters
    function GetTMVertFace(vert, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshVertUDT)      :: vert 
        integer(I8), allocatable    :: res(:)
        res = vert%face(vert%faceP(i, 1):(vert%faceP(i, 1) +  vert%faceP(i, 2) - 1))
    end function

    function GetTMCellVert(cell, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshCellUDT)      :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%vert(cell%vertP(i, 1):(cell%vertP(i, 1) + cell%vertP(i, 2) - 1))
    end function

    function GetTMCellFace(cell, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshCellUDT)      :: cell 
        integer(I8), allocatable    :: res(:)
        res = cell%face(cell%faceP(i, 1):(cell%faceP(i, 1) + cell%faceP(i, 2) - 1))
    end function

    function GetTMFaceCell(face, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshFaceUDT)      :: face 
        integer(I8), allocatable    :: res(:)
        res = face%cell(face%cellP(i, 1):(face%cellP(i, 1) + face%cellP(i, 2) - 1))
    end function

    function GetTMTubeCell(tube, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshTubeUDT)      :: tube 
        integer(I8), allocatable    :: res(:)
        res = tube%cell(tube%cellP(i, 1):(tube%cellP(i, 1) + tube%cellP(i, 2) - 1))
    end function

    function GetTMTubeFace(tube, i) result(res)
        integer(I8)                 :: i 
        class(TopomeshTubeUDT)      :: tube 
        integer(I8), allocatable    :: res(:)
        res = tube%face(tube%faceP(i, 1):(tube%faceP(i, 1) + tube%faceP(i, 2) - 1))
    end function

    function GetTMTubeBndFace(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)
        
        if (j == 1) then 
            res = tube%bndf1(tube%bndf1P(i, 1):(tube%bndf1P(i, 1) + tube%bndf1P(i, 2) - 1))
        else 
            res = tube%bndf2(tube%bndf2P(i, 1):(tube%bndf2P(i, 1) + tube%bndf2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeHighFluxBndFace(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)
        
        if (tube%hfside(i) == 1) then 
            res = tube%bndf1(tube%bndf1P(i, 1):(tube%bndf1P(i, 1) + tube%bndf1P(i, 2) - 1))
        else 
            res = tube%bndf2(tube%bndf2P(i, 1):(tube%bndf2P(i, 1) + tube%bndf2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeLowFluxBndFace(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)
        
        if (tube%hfside(i) == 1) then 
            res = tube%bndf2(tube%bndf2P(i, 1):(tube%bndf2P(i, 1) + tube%bndf2P(i, 2) - 1))
        else 
            res = tube%bndf1(tube%bndf1P(i, 1):(tube%bndf1P(i, 1) + tube%bndf1P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeHighFluxBndVert(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)

        if (tube%hfside(i) == 1) then 
            res = tube%bndv1(tube%bndv1P(i, 1):(tube%bndv1P(i, 1) + tube%bndv1P(i, 2) - 1))
        else 
            res = tube%bndv2(tube%bndv2P(i, 1):(tube%bndv2P(i, 1) + tube%bndv2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeLowFluxBndVert(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)

        if (tube%hfside(i) == 2) then 
            res = tube%bndv1(tube%bndv1P(i, 1):(tube%bndv1P(i, 1) + tube%bndv1P(i, 2) - 1))
        else 
            res = tube%bndv2(tube%bndv2P(i, 1):(tube%bndv2P(i, 1) + tube%bndv2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeBndVert(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)

        if (j == 1) then 
            res = tube%bndv1(tube%bndv1P(i, 1):(tube%bndv1P(i, 1) + tube%bndv1P(i, 2) - 1))
        else 
            res = tube%bndv2(tube%bndv2P(i, 1):(tube%bndv2P(i, 1) + tube%bndv2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeNeig(tube, i, j) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i, j
        integer(I8), allocatable    :: res(:)

        if (j == 1) then 
            res = tube%ftneig1(tube%ftneig1P(i, 1):(tube%ftneig1P(i, 1) + tube%ftneig1P(i, 2) - 1))
        else 
            res = tube%ftneig2(tube%ftneig2P(i, 1):(tube%ftneig2P(i, 1) + tube%ftneig2P(i, 2) - 1))
        end if 
    end function

    function GetTMTubeHighFluxNeig(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)

        if (tube%hfside(i) == 1) then 
            res = tube%ftneig1(tube%ftneig1P(i, 1):(tube%ftneig1P(i, 1) + tube%ftneig1P(i, 2) - 1))
        else 
            res = tube%ftneig2(tube%ftneig2P(i, 1):(tube%ftneig2P(i, 1) + tube%ftneig2P(i, 2) - 1))
        end if 

    end function

    function GetTMTubeLowFluxNeig(tube, i) result(res)
        class(TopomeshTubeUDT)      :: tube
        integer(I8), intent(in)     :: i
        integer(I8), allocatable    :: res(:)

        if (tube%hfside(i) == 2) then 
            res = tube%ftneig1(tube%ftneig1P(i, 1):(tube%ftneig1P(i, 1) + tube%ftneig1P(i, 2) - 1))
        else 
            res = tube%ftneig2(tube%ftneig2P(i, 1):(tube%ftneig2P(i, 1) + tube%ftneig2P(i, 2) - 1))
        end if 
        
    end function

    function GetTMVertFaceNeig(vert, i, tf) result(res)

        ! Description
        !============
        ! Return the 'left' and 'right' neighbouring faces for vertex tv and
        ! current face tf. If the face is not found, fn is empty. Otherwise, fn is
        ! a 1-by-2 integer array containing the face indices.

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)     :: i, tf 
        class(TopomeshVertUDT)      :: vert 
        integer(I8), allocatable    :: res(:)

        ! Auxiliary
        integer(I8)                     :: tfind
        integer(I8), allocatable        :: tvf(:), tvfe(:)

        ! Get vertex faces
        tvf = GetTMVertFace(vert, i)

        ! Find where this face is located
        tfind = findloc(tvf, tf, 1)

        ! Find the neighbours 
        if (tfind == 0) then 
            allocate(res(0))
        else
            tvfe = [tvf(size(tvf)), tvf, tvf(1)]
            res = [tvfe(tfind), tvfe(tfind+2)]
        end if 
    end function

    function AreTMFacesAdjacent(face, f1, f2) result(res)

        ! Description
        !============
        ! This routine checks if two topological mesh faces are adjacent.
        ! Faces are adjacent if they share at least one vertex.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshFaceUDT), intent(in)          :: face
        integer(I8), intent(in)                     :: f1, f2 
        logical                                     :: res 

        ! Compute
        !========
        res = .false. 
        if ((face%vert(f1, 1) == face%vert(f2, 1)) .or. &
            (face%vert(f1, 1) == face%vert(f2, 2)) .or. &
            (face%vert(f1, 2) == face%vert(f2, 1)) .or. &
            (face%vert(f1, 2) == face%vert(f2, 2))) then 
            res = .true.
        end if 

    end function

    ! ID getters
    function GetInternalFaceIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get all internal face IDs (simply faces with two cell 
        ! neighbours)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        allocate(ID(count(.not. topomesh%face%BF)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], .not. topomesh%face%BF)

    end function

    function GetBoundaryFaceIDs(topomesh) result(ID)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        allocate(ID(count(topomesh%face%BF)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], topomesh%face%BF)

    end function

    function GetSeparatrixFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are separatrix parts

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = topomesh%face%type == TMfacesepID
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetCoreFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are core parts

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = topomesh%face%type == TMfacecoreID
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetVesselFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are vessel parts (i.e. type TMfacebndID or TMfacealbndID)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = (topomesh%face%type == TMfacebndID) .or.  & 
            (topomesh%face%type == TMfacealbndID)
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)

    end function

    function GetTargetFaceIDs(topomesh) result(ID)

        ! Description
        !=============
        ! Get all face IDs that are vessel parts (i.e. type TMfacebndID)
        ! and of which at least one vertex lies on a separatrix face.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tempID, tv, tvf
        logical, allocatable, dimension(:)      :: temp

        ! Loop
        integer(I8)                 :: i, j

        ! Get
        !====
        ! Find vessel faces
        tempID = topomesh%GetVesselFaceIDs()

        ! Check which vessel faces have a separatrix vertex
        allocate(temp(size(tempID)))
        temp = .false. 
        do i = 1, size(tempID)
            ! Get face vertices
            tv = topomesh%face%vert(tempID(i), :)

            ! Check
            do j = 1, size(tv)
                tvf = topomesh%vert%GetFace(tv(j))
                if (any(topomesh%face%type(tvf) == TMfacesepID)) then 
                    temp(i) = .true.
                    exit 
                end if 
            end do 
        end do 

        ! Get indices
        allocate(ID(count(temp)))
        ID = pack(tempID, temp)

    end function

    function GetLastFluxSurfaceFaceIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get faces that form last flux surfaces (e.g. PF or SOL
        ! boundary) - these are boundary faces that are not 
        ! vessel faces

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        logical, allocatable        :: temp(:)

        ! Loop
        integer(I8)                 :: k 

        ! Get
        !====
        temp = (topomesh%face%type /= TMfacebndID) .and. &
            topomesh%face%BF 
        allocate(ID(count(temp)))
        ID = pack([(k, k = 1, topomesh%face%ntot)], temp)
        
    end function

    function GetCoreCellIDs(topomesh) result(ID)

        ! Description
        !============
        ! Get the IDs of cells that are core regions. These either have
        ! a core boundary as a face, or they have a maximum/minimum as
        ! a vertex.

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                  :: topomesh 
        integer(I8), allocatable            :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv 
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                         :: i, k

        ! Loop
        !=====
        ! Initialize
        allocate(keepind(topomesh%cell%ntot))
        keepind = .false. 

        ! Loop
        do i = 1, topomesh%cell%ntot 
            ! Check faces & vert
            tcf = topomesh%cell%GetFace(i)
            tcv = topomesh%cell%GetVert(i)

            ! Check
            if (any(topomesh%face%type(tcf) == TMfacecoreID)) then 
                keepind(i) = .true. 
            end if 
            if ((any(topomesh%vert%type(tcv) == TMvertexmaxID)) .or. &
                (any(topomesh%vert%type(tcv) == TMvertexminID))) then 
                    keepind(i) = .true.
            end if 
        end do 

        ! Extract IDs
        allocate(ID(count(keepind)))
        ID = pack([(k, k = 1, topomesh%cell%ntot)], keepind)

    end function

    function GetWideGridCellIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the cells that are consider to be a 
        ! 'wide grid' cell. This includes basically all cells that 
        ! are not next to a separatrix, or cells that have an extremum.
        ! Note that if no core boundaries are added and the grid extends
        ! to an extremum, even the typical 'core' cell is a wide grid
        ! cell!

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), allocatable, dimension(:)  :: ID 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tcf, tcv
        logical, allocatable, dimension(:)      :: keepind

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Keeper
        allocate(keepind(topomesh%cell%ntot))
        keepind = .true. 

        ! Determine
        !==========
        do i = 1, topomesh%cell%ntot 
            ! Get faces and vertices of this cell
            tcf = topomesh%cell%GetFace(i)
            tcv = topomesh%cell%GetVert(i)

            ! Any extrema?
            if ((any(topomesh%vert%type(tcv) == TMvertexmaxID)) .or. & 
                (any(topomesh%vert%type(tcv) == TMvertexminID))) then 
                keepind(i) = .false. 
                cycle
            end if 

            ! Any separatrix parts? 
            if (.not. (any(topomesh%face%type(tcf) == TMfacesepID))) then 
                keepind(i) = .false.
            end if 
        end do 

        ! Get IDs
        allocate(ID(count(keepind)))
        ID = pack([(i, i = 1, topomesh%cell%ntot)], keepind)

    end function

    function GetSeparatrixFluxSurfaceIDs(topomesh) result(fsID)

        ! Description
        !============
        ! Get separatrix flux surface IDs by checking fsIDs of all 
        ! separatrix faces

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), allocatable, dimension(:)  :: fsID
        
        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: fsIDs
        logical, allocatable, dimension(:)      :: issepface

        ! Loop

        ! Loop over faces
        !================
        issepface = topomesh%face%type == TMfacesepID
        allocate(fsIDs(count(issepface)))
        fsIDs = pack(topomesh%face%fsID, issepface)
        call Unique(fsIDs, fsID)

    end function

    function GetStrikePointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the vertex indices that are strike
        ! points. Strike points are boundary points on a connected
        ! separatrix component that contains an X-point. A disconnected
        ! boundary-to-boundary contour at a separatrix flux value is not
        ! a strike-point separatrix.
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: sepfaceIDs, &
            facevertices, sepfaces, component, stack
        logical, allocatable, dimension(:)      :: isstrikepoint, visited
        logical                                 :: hassaddle

        ! Loop
        integer(I8)                 :: i, j, k, l, currentv, ncomponent, &
            nstack

        ! Get
        !====
        ! Find separatrix faces
        sepfaceIDs = topomesh%GetSeparatrixFaceIDs()

        ! Traverse each connected separatrix component. Only boundary
        ! vertices belonging to a component with a saddle are strike
        ! points.
        allocate(isstrikepoint(topomesh%vert%ntot), &
            visited(topomesh%vert%ntot), component(topomesh%vert%ntot), &
            stack(topomesh%vert%ntot))
        isstrikepoint = .false.
        visited = .false.
        do i = 1, size(sepfaceIDs)
            facevertices = topomesh%face%vert(sepfaceIDs(i), :)
            if (visited(facevertices(1))) cycle

            nstack = 1
            stack(nstack) = facevertices(1)
            visited(facevertices(1)) = .true.
            ncomponent = 0
            hassaddle = .false.

            do while (nstack > 0)
                currentv = stack(nstack)
                nstack = nstack - 1
                ncomponent = ncomponent + 1
                component(ncomponent) = currentv
                if (topomesh%vert%type(currentv) == &
                    TMvertexsaddleID) hassaddle = .true.

                sepfaces = topomesh%vert%GetFace(currentv)
                sepfaces = pack(sepfaces, &
                    topomesh%face%type(sepfaces) == TMfacesepID)
                do k = 1, size(sepfaces)
                    facevertices = topomesh%face%vert(sepfaces(k), :)
                    do l = 1, size(facevertices)
                        if (.not. visited(facevertices(l))) then
                            nstack = nstack + 1
                            stack(nstack) = facevertices(l)
                            visited(facevertices(l)) = .true.
                        end if
                    end do
                end do
            end do

            if (hassaddle) then
                do k = 1, ncomponent
                    currentv = component(k)
                    if (topomesh%vert%type(currentv) == &
                        TMvertexbndID) isstrikepoint(currentv) = .true.
                end do
            end if
        end do

        ! Get indices
        allocate(ID(count(isstrikepoint)))
        ID = pack([(j, j = 1, topomesh%vert%ntot)], isstrikepoint)

    end function

    function GetXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the X-point (saddle point) IDs

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        logical, allocatable, dimension(:)      :: isXP

        ! Loop
        integer(I8)                             :: k

        ! Find X-points
        !==============
        isXP = topomesh%vert%type == TMvertexsaddleID
        allocate(ID(count(isXP)))
        ID = pack([(k, k = 1, size(isXP))], isXP)

    end function

    function GetOPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns all the O-points (extrema) point) IDs

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        logical, allocatable, dimension(:)      :: isOP

        ! Loop
        integer(I8)                             :: k

        ! Find X-points
        !==============
        isOP = (topomesh%vert%type == TMvertexminID) .or. &
            (topomesh%vert%type == TMvertexmaxID)
        allocate(ID(count(isOP)))
        ID = pack([(k, k = 1, size(isOP))], isOP)

    end function

    function GetPrimaryXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the vertex IDs of all x-points that 
        ! lie close to a core (close in the sense that there is no other
        ! x-point between this x-point and the next o-point). This can 
        ! only occur if the x-point is present in a core region's 
        ! vertices. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        integer(I8), allocatable, dimension(:)  :: ID

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: corecells, tv

        ! Loop
        integer(I8)                             :: i

        ! Initialize
        !===========
        ! Allocate
        allocate(ID(0))

        ! Unpack
        associate(&
            cell        => topomesh%cell,   &
            vert        => topomesh%vert    &
            )

        ! Determine primary x-points
        !===========================
        ! Determine core cells
        corecells = topomesh%GetCoreCellIDs()
        do i = 1, size(corecells)
            ! Get vertices
            tv = cell%GetVert(corecells(i))

            ! Check which ones are x-points (if any) and add
            ID = [ID, pack(tv, vert%type(tv) == TMvertexsaddleID)]
        end do 

        ! Housekeeping
        !=============
        end associate

    end function

    function GetStrikePointXPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the x-point ID for each strike point 
        ! as returned by GetStrikePointIDs (a strike point can only 
        ! have one X-point normally speaking). To determine this, we 
        ! start from each strike point and keep tracing the aligned 
        ! field line it belongs to until an x-point is reached (normally, 
        ! this is already reached in the first face but ok)
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8)                             :: nsp, &
            nextv, tsp
        integer(I8), allocatable, dimension(:)  ::sp, tspf, tfvert
        logical, allocatable, dimension(:)      :: isfacefound

        ! Loop
        integer(I8)                 :: i

        ! Initialize
        !===========
        ! Get strike points
        sp = topomesh%GetStrikePointIDs()
        nsp = size(sp)
        
        ! Initialize array
        allocate(ID(nsp), isfacefound(topomesh%face%ntot))
        ID = 0

        ! Associate
        associate(&
            vert    => topomesh%vert,   &
            face    => topomesh%face    &
            )

        ! Loop
        !=====
        do i = 1, nsp
            ! Initialize
            isfacefound = .false. 

            ! Unpack 
            tsp = sp(i)

            ! Loop until vertex is found or until we reach the next 
            ! boundary again
            do while (.true.)

                ! Get faces of this point
                tspf = vert%GetFace(tsp)

                ! Get separatrix faces that have not yet been found
                tspf = pack(tspf, (face%type(tspf) == TMfacesepID) .and. &
                    (.not. isfacefound(tspf)))
                
                ! Sanity checks
                if (size(tspf) < 1) then 
                    call gdErrorHandler('GetStrikePointXPointIDs: could not ' // & 
                        'find next separatrix segment, unexpected')
                end if 
                if (size(tspf) > 1) then 
                    call gdErrorHandler('GetStrikePointXPointIDs: found more ' // & 
                        'than one separatrix segment, unexpected')
                end if 

                ! Set separatrix face as found
                isfacefound(tspf) = .true. 

                ! Get the next vertex
                tfvert = face%vert(tspf(1), :)
                if (tfvert(1) == tsp) then 
                    nextv = tfvert(2)
                elseif (tfvert(2) == tsp) then 
                    nextv = tfvert(1)
                else
                    call gdErrorHandler('GetStrikePointXPointIDs: something ' // & 
                        'wrong with topomesh interconnections, could not ' // & 
                        'find current vertex in face vertices')
                end if 

                ! Check the point type
                if (vert%type(nextv) == TMvertexsaddleID) then 
                    ! Found, add and exit
                    ID(i) = nextv 
                    exit 
                elseif (vert%type(nextv) == TMvertexbndID) then 
                    ! Shouldn't happen, error
                    call gdErrorHandler('GetStrikePointXPointIDs: could not ' // & 
                        'find X-point of strike point')
                else
                    ! Need to look further along the separatrix. Set the
                    ! current vertex as the next vertex
                    tsp = nextv 
                end if 
            end do 
        end do

        ! Housekeeping
        !=============
        end associate

    end function

    function GetClosedContourTangencyPointIDs(topomesh) result(ID)

        ! Description
        !============
        ! This function returns the tangency point IDs (type 2) that 
        ! have a contour which closes upon itself. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)          :: topomesh 
        integer(I8), allocatable    :: ID(:)

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tv, &
            sortind, tf
        integer(I8), allocatable, dimension(:, :)   :: tfvert
        logical, allocatable, dimension(:)      :: isbranchingpolygon, &
            ispolygonstart

        ! Loop
        integer(I8)                 :: i, k

        ! Initialize
        !===========
        ! Allocate
        allocate(ID(0))

        ! Associate
        associate(&
            face        => topomesh%face,   &
            vert        => topomesh%vert    &
            )

        ! Loop
        !=====
        do i = 1, vert%ntot
            if (vert%type(i) == TMvertextp2ID) then 
                ! Get all topomesh faces with this flux surface ID
                allocate(tf(count(face%fsID == vert%fsID(i))))
                tf = pack([(k, k = 1, face%ntot)], face%fsID == vert%fsID(i))

                ! Check by sorting vertices 
                tfvert = face%vert(tf, :)
                allocate(sortind(size(tf)), ispolygonstart(size(tf)), &
                    isbranchingpolygon(size(tf)))
                call SortPolygonEdges(tfvert, size(tf), sortind, ispolygonstart, &
                    isbranchingpolygon)
                tfvert = tfvert(sortind, :)

                ! Sanity checks
                if (count(ispolygonstart) == 0) then 
                    ! Possible if aligned vessel parts were inserted, 
                    ! skip
                    deallocate(tf, sortind, ispolygonstart, isbranchingpolygon)
                    cycle
                elseif (count(ispolygonstart) > 1) then 
                    ! Possible if aligned vessel parts were inserted, 
                    ! skip
                    deallocate(tf, sortind, ispolygonstart, isbranchingpolygon)
                    cycle
                end if 
                if (count(isbranchingpolygon) /= 0) then 
                    call gdErrorHandler('GetClosedContourTangencyPointIDs: ' // &
                        'found branching polygons, unexpected')
                end if 

                ! Check first and last edge
                call SetDiff(tfvert(1, :), tfvert(size(tf), :), tv)
                if (size(tv) == 0) then 
                    ! Found, add
                    ID = [ID, i]
                end if 

                ! Housekeeping
                deallocate(tf, sortind, ispolygonstart, isbranchingpolygon)

            end if 
        end do

        ! Housekeeping
        !=============
        end associate

    end function
 
    ! Topological mesh identification 
    function IdentifyTopologicalMeshType(topomesh) result(TMlabel)

        ! Description
        !============
        ! This function attempts to identify the topological mesh type
        ! and returns the topological mesh identification number 
        ! (see mod_definitions of the definition thereof). This is not 
        ! really required by goat itself, but is rather a convenience
        ! tool for later simulation etc and face label mapping. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)              :: topomesh
        integer(I8)                     :: TMlabel

        ! Auxiliary
        logical                                 :: issinglenull, &
            isdoublenull, islinear, islimiter
        integer(I8)                             :: nxp, nop, nwgc, nsp, &
            ncc, ntp
        integer(I8), allocatable, dimension(:)  :: xp, op, wgc, sp, cc, tp

        ! Initialize
        !===========
        ! Set initial label to default label
        TMlabel = TMTopGeneral

        ! Get all X-points
        xp = topomesh%GetXPointIDs()
        nxp = size(xp)

        ! Get all O-points
        op = topomesh%GetOPointIDs()
        nop = size(op)

        ! Strike points
        sp = topomesh%GetStrikePointIDs()
        nsp = size(sp)

        ! Closed contour tangency points
        tp = topomesh%GetClosedContourTangencyPointIDs()
        ntp = size(tp)

        ! Get all wide grid cells
        wgc = topomesh%GetWideGridCellIDs()        
        nwgc = size(wgc)

        ! Core cells
        cc = topomesh%GetCoreCellIDs()
        ncc = size(cc)

        ! Initialize
        issinglenull = .true.
        isdoublenull = .true.
        islinear     = .true. 
        islimiter    = .true.

        ! Linear case
        !============
        ! X-point and O-point checks
        if (nxp /= 0 .or. nop /= 0 .or. ntp /= 0 .or. ncc /= 0) then
            islinear = .false. 
        end if 

        ! Limiter case
        !=============
        ! A limiter mesh has no X/O points, but it does have closed-flux
        ! core topology. In practice this appears either as closed-contour
        ! tangency points, or as core cells bounded by the core flux surface.
        if (nxp /= 0 .or. nop /= 0) then
            islimiter = .false.
        end if
        if (ntp == 0 .and. ncc == 0) then
            islimiter = .false.
        end if

        ! Single null
        !============        
        ! X-point and O-point checks
        if (nxp /= 1 .or. nop /= 0) then 
            issinglenull = .false.
        end if 

        ! Strike point checks
        if (nsp /= 2) then 
            issinglenull = .false. 
        end if 

        ! Number of wide and narrow grid cells
        !if (nwgc > 0) then 
        !    issinglenull = .false.
        !elseif (topomesh%cell%ntot /= 3) then
        !    issinglenull = .false.
        !end if 

        ! Core cells
        !if (ncc /= 1) then 
        !    issinglenull = .false.
        !end if 

        ! Double null
        !============
        ! X-point and O-point checks
        if (nxp /= 2 .or. nop /= 0) then 
            isdoublenull = .false.
        end if 

        ! Strike point checks
        !if (nsp /= 4) then 
        !    isdoublenull = .false. 
        !end if 

        ! Number of wide and narrow grid cells
        !if (nwgc > 0) then 
        !    isdoublenull = .false.
        !elseif (topomesh%cell%ntot /= 3) then
        !    isdoublenull = .false.
        !end if 

        ! Core cells
        !if (ncc /= 2) then 
        !    isdoublenull = .false.
        !end if 

        ! Determine flag
        !===============
        ! Sanity checks
        if (count([issinglenull, isdoublenull, islinear, islimiter]) > 1) then
            ! Probably we missed something in the definition then
            print *, 'IdentifyTopologicalMeshType: multiple topologies ' // & 
                'appear valid, this is likely a bug. Setting flag to ' // & 
                'general flag...'
        elseif (islinear) then 
            TMlabel = TMTopL
        elseif (islimiter) then
            TMlabel = TMTopLM
        elseif (issinglenull) then 
            TMlabel = TMTopSN
        elseif (isdoublenull) then 
            TMlabel = TMTopDN
        end if 

    end function

    ! Topological metric computations
    subroutine GetTMFacePsiValueDistribution(topomesh, fieldtracer, &
        faceID, psi, psi_orig)

        ! Description
        !============
        ! This function returns the monotonized psi value distribution 
        ! on a topological mesh face. If it's an aligned face, then 
        ! the distribution is uniform and equal to the flux surface 
        ! ID of the face. If it is not an aligned face, then the 
        ! distribution is evaluated using the fieldtracer at all 
        ! face points and made monotonically increasing or decreasing, 
        ! depending on the values in the face vertices. 

        ! Note: for aligned boundary faces, this routine will give 
        ! inaccurate results, since there is likely some variation along
        ! the face. Thread carefully in that case. 

        ! Note: we also return the original, non-monotonized distribution 
        ! in psi_orig

        ! Declare variables
        !==================
        ! Arguments
        type(TopomeshUDT), intent(in)                   :: topomesh
        class(ContourtracerUDT), intent(in)             :: fieldtracer
        integer(I8), intent(in)                         :: faceID
        real(R8), allocatable, dimension(:), intent(out)    :: psi, psi_orig

        ! Auxiliary
        real(R8), allocatable, dimension(:)             :: xf, yf

        ! Loop
        integer(I8)                                     :: i

        ! Compute
        !========
        ! Get face coordinates
        xf = topomesh%face%x(faceID)%Get()
        yf = topomesh%face%y(faceID)%Get()

        ! Initialize
        allocate(psi(size(xf)))
        psi = 0.0_R8
        
        ! Hedge for aligned faces - unless they're aligned boundary faces
        if ((topomesh%face%fsID(faceID) /= 0) .and. (topomesh%face%type(faceID) /= TMfacealbndID)) then 
            ! Set uniform distribution and return
            psi = topomesh%fsfval%Get(topomesh%Face%fsID(faceID))
            psi_orig = psi
            return 
        end if

        ! If we got here, it is a non-aligned face
        psi = fieldtracer%Evaluate(xf, yf)
        psi_orig = psi
        psi(1) = topomesh%vert%fval(topomesh%face%vert(faceID, 1))
        psi(size(psi)) = topomesh%vert%fval(topomesh%face%vert(faceID, 2))

        ! Monotonize
        if (psi(1) < psi(size(psi))) then 
            ! Increasing psi
            do i = 2, size(psi)
                if (psi(i) < psi(i-1)) then 
                    psi(i) = psi(i-1)
                end if 
            end do  
        else
            ! Decreasing psi
            do i = 2, size(psi)
                if (psi(i) > psi(i-1)) then 
                    psi(i) = psi(i-1)
                end if 
            end do
        end if


    end subroutine 

    function GetTMFaceRadialLengthDistribution(topomesh, fieldtracer, &
        magneticField, faceID) result(dlcrad)

        ! Description
        !============
        ! This function computes the (absolute) length distribution of a face along
        ! the radial direction. This is determined as the sum of the 
        ! radial lengths of the face's edges. 

        ! Note: if parts of the face exhibit non-monotonous behavior 
        ! in terms of psi value, the length of these parts is set to zero

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)      :: topomesh 
        class(ContourtracerUDT), intent(in) :: fieldtracer
        integer(I8), intent(in)             :: faceID
        type(MagneticFieldUDT), intent(in)  :: magneticField 
        real(R8), allocatable, dimension(:) :: dlcrad 

        ! Auxiliary
        real(R8), allocatable, dimension(:) :: x, y, dx, dy, bx, by, &
            xf, yf, bn, psi

        ! Loop
        integer(I8)                         :: i

        ! Compute
        !========
        ! Get face coordinates
        x = topomesh%face%x(faceID)%Get()
        y = topomesh%face%y(faceID)%Get()
        
        ! Initialize
        allocate(dlcrad(size(x)))
        dlcrad = 0.0_R8

        ! Compute edge center coordinates and lengths, and psi gradient
        dx = x(2:) - x(1:size(x)-1)
        dy = y(2:) - y(1:size(y)-1)
        xf = 0.5*(x(2:) + x(1:size(x)-1))
        yf = 0.5*(y(2:) + y(1:size(y)-1))
        allocate(bx(size(xf)), by(size(xf)), psi(size(x)))
        call magneticField%interp%Evaluate(x , y , 0, 0, psi)
        call magneticField%interp%Evaluate(xf, yf, 1, 0, bx)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, by)
        bn = sqrt(bx**2 + by**2)
        bx = bx/bn
        by = by/bn

        ! Monotonize
        if (psi(1) < psi(size(psi))) then 
            ! Increasing psi
            do i = 2, size(psi)
                if (psi(i) < psi(i-1)) then 
                    psi(i) = psi(i-1) ! also keep psi the same!
                    dlcrad(i) = dlcrad(i-1)
                else
                    dlcrad(i) = dlcrad(i-1) + abs(bx(i-1)*dx(i-1) + by(i-1)*dy(i-1))
                end if 
            end do  
        else
            ! Decreasing psi
            do i = 2, size(psi)
                if (psi(i) > psi(i-1)) then 
                    psi(i) = psi(i-1)
                    dlcrad(i) = dlcrad(i-1)
                else
                    dlcrad(i) = dlcrad(i-1) + abs(bx(i-1)*dx(i-1) + by(i-1)*dy(i-1))
                end if 
            end do
        end if

        ! Housekeeping
        deallocate(bx, by, psi)

    end function

    subroutine GetTMTubePsiLimitsTA(tmadaptor, topomesh, tubeID, &
        psimin, psimax, includealbndin)

        ! Description
        !============
        ! Get the tube psi value limits (i.e. the maximal value of 
        ! psi at the low psi value bound and the minimal value of psi 
        ! at the high value bound - then all contours traced in between
        ! those bounds should lie nicely within the tube). If, for 
        ! whatever reason, the psi values overlap, then psimin will be 
        ! larger than psimax and a message will be shown. 

        ! Note: contributions of aligned boundary faces are no longer
        ! included here by default. The reason is that this routine is used to
        ! compute merging criteria, and aligned boundaries cannot be 
        ! merged away. We should find a more elegant solution at some 
        ! point... If aligned boundaries should be included, one can 
        ! use the optional argument includealbndin and set it to .true. 

        ! Note: this routine is very similar to the (old) GetTMTubesPsiLimitsTA
        ! routine, with the exception that the psi values are now evaluated
        ! on the true vertex locations of each boundary. This is necessary
        ! to avoid wrong limits (and hence possible issues downstream)
        ! because tangency points etc are removed by topomesh cleaning 
        ! routines. It also provides a more discretely consistent 
        ! approach. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT), intent(in)   :: tmadaptor
        type(TopomeshUDT), intent(in)           :: topomesh 
        integer(I8), intent(in)                 :: tubeID
        real(R8), intent(out)                   :: psimin, psimax
        logical, intent(in), optional           :: includealbndin

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tf1, tf2, tv1, tv2
        real(R8), allocatable, dimension(:)     :: psi1, psi2
        logical                                 :: includealbnd

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Check if we should include aligned boundaries 
        if (present(includealbndin)) then 
            includealbnd = includealbndin
        else
            includealbnd = .false. 
        end if 

        ! Unpack for ease
        associate(&
            tube        => topomesh%tube,   &
            vert        => topomesh%vert,   &
            face        => topomesh%face    &
            )

        ! Compute bounds
        !===============
        ! Get faces and vertices at both sides
        tf1 = tube%GetBndFace(tubeID, 1)
        tf2 = tube%GetBndFace(tubeID, 2)
        tv1 = tube%GetBndVert(tubeID, 1)
        tv2 = tube%GetBndVert(tubeID, 2)

        if (includealbnd) then 
            ! Keep only vertices and faces with non-zero ID and any 
            ! vertices that are not type 1 tangency points
            tf1 = pack(tf1, (face%fsID(tf1) /= 0))
            tf2 = pack(tf2, (face%fsID(tf2) /= 0))
            tv1 = pack(tv1, vert%fsID(tv1) /= 0) ! keep all vertices with flux surface ID - flux surface value of aligned boundaries may be inaccurate
            tv2 = pack(tv2, vert%fsID(tv2) /= 0)
        else
            ! Keep only vertices and faces with non-zero ID and remove any
            ! faces that are aligned boundaries and any vertices that are not 
            ! type 1 tangency points
            tf1 = pack(tf1, (face%fsID(tf1) /= 0) .and. (face%type(tf1) /= TMfacealbndID))
            tf2 = pack(tf2, (face%fsID(tf2) /= 0) .and. (face%type(tf2) /= TMfacealbndID))
            tv1 = pack(tv1, vert%type(tv1) == TMvertextp1ID .or. &
                vert%type(tv1) == TMvertexminID .or. vert%type(tv1) == TMvertexmaxID) !only keep type 1 TPs, max and min
            tv2 = pack(tv2, vert%type(tv2) == TMvertextp1ID .or. &
                vert%type(tv2) == TMvertexminID .or. vert%type(tv2) == TMvertexmaxID)

            ! If there are sides with only aligned boundaries, take these and
            ! issue warning
            if ((size(tf1) + size(tv1)) == 0) then 
                ! May happen in some cases
                print *, ('warning: GetTMTubePsiLimitsTA: tube has only aligned ' // &
                    'boundary faces as neighbour at side 1, including these to determine ')
                tf1 = tube%GetBndFace(tubeID, 1)
                tv1 = tube%GetBndVert(tubeID, 1)
                tf1 = pack(tf1, (face%fsID(tf1) /= 0))
                tv1 = pack(tv1, vert%fsID(tv1) /= 0)
            end if 
            if ((size(tf2) + size(tv2)) == 0) then 
                ! May happen in some cases
                print *, ('warning: GetTMTubePsiLimitsTA: tube has only aligned ' // &
                'boundary faces as neighbour at side 2, including these to determine ')
                tf2 = tube%GetBndFace(tubeID, 2)
                tv2 = tube%GetBndVert(tubeID, 2)
                tf2 = pack(tf2, (face%fsID(tf2) /= 0))
                tv2 = pack(tv2, vert%fsID(tv2) /= 0)
            end if 
        end if 

        ! Check
        if ((size(tf1) + size(tv1)) == 0) then 
            ! Should not happen
            call WriteTopologicalMesh(topomesh, 'topomesh_error')
            call gdErrorHandler('GetTMTubePsiLimitsTA: tube has only aligned ' // &
                'boundary faces as neighbour at side 1, unexpected')
        end if 
        if ((size(tf2) + size(tv2)) == 0) then 
            ! Should not happen
            call WriteTopologicalMesh(topomesh, 'topomesh_error')
            call gdErrorHandler('GetTMTubePsiLimitsTA: tube has only aligned ' // &
                'boundary faces as neighbour at side 2, unexpected')
        end if 

        ! Get psi values
        allocate(psi1(0), psi2(0))
        do i = 1, size(tv1)
            psi1 = [psi1, topomesh%fsfval%Get(vert%fsID(tv1(i)))]
        end do 
        do i = 1, size(tf1)
            psi1 = [psi1, tmadaptor%facepsi(tf1(i))%Get()]
        end do 
        do i = 1, size(tv2)
            psi2 = [psi2, topomesh%fsfval%Get(vert%fsID(tv2(i)))]
        end do 
        do i = 1, size(tf2)
            psi2 = [psi2, tmadaptor%facepsi(tf2(i))%Get()]
        end do 

        ! Check which side is low and which is high
        if (all(minval(psi1) > psi2)) then 
            ! First side is high flux boundary
            psimax = minval(psi1)
            psimin = maxval(psi2)
        elseif (all(minval(psi2) > psi1)) then 
            ! Second side is high flux boundary
            psimax = minval(psi2)
            psimin = maxval(psi1)
        else 
            print *, 'GetTMTubePsiLimitsTA: psi values seem to overlap, ' // &
                'could not determine high and low psi side. psimin will ' // &
                'be larger than psimax...'
            psimax = minval(psi1)
            psimin = maxval(psi2)
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    subroutine GetTMTubeRadialWidthTA(tmadaptor, topomesh, tubeID, &
            lrad, faceID, dlcradface, includealbndin)

        ! Description
        !============
        ! This routine returns the radial width of a topological mesh 
        ! tube, taking into account also possible differences in psi 
        ! value of the bounding flux surfaces. The radial length is 
        ! computed as follows: 
        ! - for each face of the tube, the radial length is computed
        ! - the minimal value of this radial length is determined and
        ! returned
        ! - the radial length is computed only between the psi value 
        ! bounds of the tube. The field tracer is used to determine the
        ! psi values on each face to ensure discretely consistent results
        ! - if the psi value bounds overlap, then the radial width is zero
        ! Optionally, the routine also returns the faceID with the 
        ! minimal length and its radial length distribution on all of 
        ! its vertices. 

        ! Note: we parse the option includealbndin to optionally include
        ! aligned boundary psi values when computing the length. This 
        ! may lead to zero width (we ceil the length to zero to be 
        ! nonnegative), but may be useful in some cases. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshAdaptorUDT)               :: tmadaptor
        type(TopomeshUDT), intent(in)           :: topomesh
        integer(I8), intent(in)                 :: tubeID 
        real(R8), intent(out)                   :: lrad
        integer(I8), intent(out), optional      :: faceID
        real(R8), allocatable, dimension(:), intent(out), optional  :: dlcradface
        logical, optional, intent(in)           :: includealbndin 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: tubef
        real(R8)                                :: psimin, psimax, &
            thislrad
        real(R8), allocatable, dimension(:)     :: xf, yf, psif, dlcradf, &
            lradminmax, dlradf
        logical                                 :: includealbnd

        ! Loop
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Check optional inputs
        if (present(includealbndin)) then 
            includealbnd = includealbndin
        else
            includealbnd = .false.
        end if 

        ! Initialize
        lrad = posinfval_R8()

        ! Unpack
        associate(&
            tube    => topomesh%tube,   &
            face    => topomesh%face    &
            )

        ! Compute
        !========
        ! Psi bounds
        call GetTMTubePsiLimitsTA(tmadaptor, topomesh, tubeID, &
            psimin, psimax, includealbnd)

        ! Get tube faces
        tubef = tube%GetFace(tubeID)

        ! Check
        if (psimin >= psimax) then 
            print *, 'GetTMTubeRadialWidthTA: psimin >= psimax, returning ' // &
                'zero radial width'
            lrad = 0.0_R8 
            if (present(faceID)) then 
                faceID = tubef(1)
            end if 
            if (present(dlcradface)) then 
                xf = face%x(tubef(1))%Get()
                allocate(dlcradface(size(xf)))
                dlcradface = 0.0_R8
            end if 
            return 
        end if 

        ! Loop over faces
        do i = 1, size(tubef)
            ! Get face coordinates and monotonized psi distribution 
            xf = face%x(tubef(i))%Get()
            yf = face%y(tubef(i))%Get()
            psif = tmadaptor%facepsi_mono(tubef(i))%Get()

            ! Hedge for psimin/psimax not lying on this face
            if (all(psimin > psif) .or. all(psimin < psif)) then 
                psimin = minval(psif)
            end if 
            if (all(psimax > psif) .or. all(psimax < psif)) then 
                psimax = maxval(psif)
            end if 
            !psimin = max(psimin, minval(psif))
            !psimax = min(psimax, maxval(psif))

            ! Compute the monotonized face radial length distribution
            dlcradf = tmadaptor%facedlcrad(tubef(i))%Get()

            ! Compute radial length coordinate where psimin and psimax occur
            if (psif(1) > psif(size(psif))) then 
                ! switch for interpolation
                psif = psif(size(psif):1:-1)
                dlcradf = dlcradf(size(dlcradf):1:-1)
                call Interpolate1D([psimin, psimax], lradminmax, psif, dlcradf)

                ! Switch back
                dlcradf = dlcradf(size(dlcradf):1:-1)
                psif = psif(size(psif):1:-1)
            else
                call Interpolate1D([psimin, psimax], lradminmax, psif, dlcradf)
            end if 
            dlradf = dlcradf(2:) - dlcradf(1:size(dlcradf)-1)

            ! Check
            if (any(isnan(lradminmax))) then 
                call gdErrorHandler('GetTMTubeRadialWidth: nans detected ' // & 
                    'when evaluating radial length at psi boundaries of ' // &
                    'tube. This indicates these psi values are not present ' // & 
                    'on the face, unexpected')
            end if 
            
            ! Adjust dlcrad
            dlcradf = 0.0_R8 ! Rebuild
            where (psif < psimin) psif = psimin
            where (psif > psimax) psif = psimax 
            if (psif(1) < psif(size(psif))) then 
                ! Psi is ascending, monotonize
                do j = 2, size(psif)
                    if ((psif(j) <= psif(j-1))) then 
                        psif(j) = psif(j-1)
                        dlcradf(j) = dlcradf(j-1)
                    else
                        dlcradf(j) = dlcradf(j-1) + dlradf(j-1)
                    end if 
                end do  
            else
                ! Psi is descending, monotonize
                do j = 2, size(psif)
                    if ((psif(j) >= psif(j-1))) then 
                        psif(j) = psif(j-1)
                        dlcradf(j) = dlcradf(j-1)
                    else
                        dlcradf(j) = dlcradf(j-1) + dlradf(j-1)
                    end if 
                end do  
            end if

            ! Compute radial length
            thislrad = dlcradf(size(dlcradf))

            ! Check
            if (thislrad < lrad) then 
                lrad = thislrad
                if (present(faceID)) then 
                    faceID = tubef(i)
                end if 
                if (present(dlcradface)) then 
                    dlcradface = dlcradf
                end if 
            end if 
        end do 

        ! Housekeeping
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                 TOPOLOGICAL MESH CELL OPERATORS                  !
    !------------------------------------------------------------------!
    ! Cell face sorter (CW or CCW)
    subroutine SortTopologicalMeshCellFaces(topomesh, cellind, &
        sortedfaces, sortedvertices, doflip)

        ! Description
        !============
        ! This routine sorts the faces of a cell in arbitrary direction
        ! (clockwise or counter-clockwise) and returns the sorted face
        ! IDs (sortedfaces) and an indication whether the face 
        ! coordinates should be flipped when querying the face from 
        ! the topological mesh. This may be useful to post-process 
        ! cells etc. Additionally, the vertices of the cell are also
        ! returned (closed form, so first and last are the same).

        ! Note: by construction, the cell faces should already be 
        ! sorted correctly! 

        ! Note: we should hedge for possible duplicate faces since 
        ! disc regions are allowed

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: cellind 
        integer(I8), allocatable, intent(out)   :: sortedfaces(:), &
            sortedvertices(:)
        logical(I8), allocatable, intent(out)   :: doflip(:)

        ! Auxiliary
        integer(I8)                             :: startvert, sf
        integer(I8), allocatable, dimension(:)  :: tf, tfv, tdiff

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Associate
        associate(&
            cell        => topomesh%cell,   &
            face        => topomesh%face)

        ! Unpack 
        sortedfaces = GetTMCellFace(cell, cellind)
        sortedvertices = GetTMCellVert(cell, cellind)
        tf = sortedfaces

        ! Initialize
        if (allocated(doflip)) then 
            deallocate(doflip)
        end if 
        allocate(doflip(size(sortedfaces)))
        doflip = .false.

        ! Determine starting face as non-duplicate face (then all the rest
        ! should go fine in principle)
        tdiff = tf(2:size(tf)) - tf(1:size(tf)-1)
        sf = findloc(tdiff /= 0, .true., 1)
        if (sf == 1) then 
            ! Nothing to see here, move along
        elseif (sf == 2) then 
            ! Shift the first face backwards
            tf = [tf(2:size(tf)), tf(1)]
        else
            ! Unexpected, throw error
            call gdErrorHandler('SortTopologicalMeshCellFaces: more than two faces ' // & 
                'appear to be the same for this cell, unexpected')
        end if 

        ! Loop
        !=====
        if (any(face%vert(tf(1), 1) == face%vert(tf(2), :))) then 
            startvert = face%vert(tf(1), 1)
            doflip(1) = .true. 
        elseif (any(face%vert(tf(1), 2) == face%vert(tf(2), :))) then 
            startvert = face%vert(tf(1), 2)
            doflip(1) = .false. 
        else 
            call gdErrorHandler('SortTopologicalMeshCellFaces: unexpected error')
        end if

        do i = 2, size(tf)
            tfv = face%vert(tf(i), :)
            if (tfv(2) == startvert) then 
                doflip(i) = .true. 
                startvert = tfv(1)
            elseif (tfv(1) == startvert) then 
                doflip(i) = .false.
                startvert = tfv(2)
            else 
                call gdErrorHandler('SortTopologicalMeshCellFaces: unexpected error')
            end if 
        end do 
        


        ! Housekeeping
        !=============
        end associate 

    end subroutine

    ! Full cell polygon coordinate constructor
    subroutine ConstructTopologicalMeshCellPolygon(topomesh, cellind, &
        xcda, ycda)

        ! Description
        !============
        ! Construct the cell polygon by appending the face polygon
        ! coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), intent(in)                 :: cellind
        type(RealDynamicArrayUDT), intent(out)  :: xcda, ycda

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: sortedfaces, &
            sortedvertices
        real(R8), allocatable, dimension(:)     :: tempx, tempy 
        logical, allocatable, dimension(:)      :: doflip

        ! Loop
        integer(I8)                             :: i 

        ! Initialize
        !===========
        ! Initialize dynamic arrays
        xcda = ConstructRealDynamicArray()
        ycda = ConstructRealDynamicArray()

        ! Construct polygon
        !==================
        ! Sort the cell faces
        call SortTopologicalMeshCellFaces(topomesh, cellind, sortedfaces, &
            sortedvertices, doflip)

        ! Loop over the faces 
        do i = 1, size(sortedfaces)
            ! Get face coordinates
            tempx = topomesh%face%x(sortedfaces(i))%Get()
            tempy = topomesh%face%y(sortedfaces(i))%Get()

            ! Check
            if (doflip(i)) then 
                call xcda%Append(tempx(size(tempx):2:-1)) ! avoid double coordinates
                call ycda%Append(tempy(size(tempy):2:-1))
            else
                call xcda%Append(tempx(2:size(tempx))) ! avoid double coordinates
                call ycda%Append(tempy(2:size(tempy)))
            end if 
        end do 

    end subroutine 

    ! Closed polygon constructor from topomesh cells (based on cell 
    ! addition algorithm)
    function GetClosedPolygonFromTopomeshFaces(topomesh, faceID) &
        result(sortedfaces)

        ! Description
        !============
        ! This routine attempts to form a closed polygon by sorting the 
        ! given face IDs in 'faceID' based on the vertices and the faces
        ! they connect to in the topomesh. The main advantage of this
        ! routine is that there is branching polygon support, i.e. we 
        ! do a similar approach as when constructing grid cells, relying
        ! on how faces have been sorted per vertex. If no single closed
        ! polygon that contains all given faces could be formed, 
        ! the result is an empty array. Otherwise, it is the sorted 
        ! array of faceID with subsequent faces. Standard assumptions on 
        ! faces hold (e.g. no faces with twice the same vertices, unique
        ! vertex pairs per face, ...)

        ! Algorithm
        !==========
        ! We apply a similar algorithm as used for general cell 
        ! construction, but now we only allow to use faces that have 
        ! been given. Furthermore, we trace from the same starting face
        ! in both directions. 

        ! 0)    Take a starting face from the given face IDs and set the 
        !       current turning direction. Take a starting vertex.  
        ! 1)    Take the next vertex of the current face
        ! 2)    If the vertex has only one other face that is available, 
        !       take that one and go to 1) if it is not the same as the 
        !       start face. If it is the startface, then exit. If the 
        !       vertex has multiple faces, then see if the turning direction
        !       face is available. If that's the case, take the face and go to 
        !       1). If that's not the case, exit (no cell could be found
        !       for this turning direction) and go to 3)
        ! 3)    If only one direction was tried, restart from 0 with the
        !       same face but now doing the other turning direction. If
        !       both turning directions have been considered, exit. 

        ! Declare variables
        !==================
        ! Arguments 
        class(TopomeshUDT)                      :: topomesh 
        integer(I8), dimension(:), intent(in)   :: faceID
        integer(I8), allocatable, dimension(:)  :: sortedfaces

        ! Auxiliary 
        integer(I8)                             :: tf, startface, & 
            turndirection, tfv(1:2), tvind, tv, startvert, &
            starttvind, nf, nfv(1:2)
        integer(I8), allocatable, dimension(:)  :: fc, disccellvert, &
            nfvfn, tcf, tcv, faceneig1, faceneig2, faceneig, tempf, diffc, &
            vertface

        logical                                 :: istfv(1:2)
        logical, allocatable                    :: hasturned1(:, :), &
            hasturned2(:, :), donotstartfromface(:), considerface(:)

        type(IntegerDynamicArrayUDT), allocatable   :: cellvert(:), &
            cellface(:)
        type(IntegerDynamicArrayUDT)                :: thiscellvert, &
            thiscellface

        ! Loop 
        integer(I8)                             :: i, cc 

        ! Initialize
        !===========
        ! Output
        allocate(sortedfaces(0))

        ! Checks
        if (size(faceID) == 0) then 
            return 
        end if 

        ! Unpack
        associate( & 
            vert        => topomesh%vert,   &
            face        => topomesh%face,   &
            cell        => topomesh%cell    &
            )

        ! Allocate
        allocate(fc(face%ntot), hasturned1(face%ntot, 2), &
            hasturned2(face%ntot, 2), donotstartfromface(face%ntot), &
            cellvert(0), cellface(0), considerface(face%ntot))

        ! Check which faces to consider
        considerface = .false.
        considerface(faceID) = .true. 

        ! Initialize face counters
        fc = 2 ! Note: here boundary faces may be passed by twice without problems
        where (.not. considerface) fc = 0
            
        ! Initialize turn checkers
        hasturned1 = .false.
        hasturned2 = .false.

        ! Initialize cell counter
        cc = 0

        ! Check if there are any disc-type cells
        disccellvert = findloc(vert%faceP(:, 2), 1_I8)

        ! Check if there are any faces with only one adjacent face on each side -
        ! these faces shouldn't be started from, as one cannot determine the
        ! turning direction (and several other faces should remain that can be
        ! started from)
        donotstartfromface = .false. 
        do i = 1, face%ntot
            nfvfn = vert%faceP(face%vert(i, :), 2)
            if (all(nfvfn == 2)) then 
                donotstartfromface(i) = .true.
            end if 
        end do 

        ! Loop
        !=====
        turndirection = 0
        do while (.true.) ! we basically loop over the turndirection though...

            ! Turning direction (0: not found, 1: first neighbour, 2: second
            ! neighbour)
            turndirection = turndirection + 1
            if (turndirection > 2) then 
                exit 
            end if 
            
            ! Initialize cell faces & vertices
            allocate(tcv(0), tcf(0))

            ! Find the next face (any next internal face)
            tf = findloc( (fc > 0) .and. (.not. face%BF) .and. &
                (.not. donotstartfromface), .true., 1)
            
            ! Check
            if (tf == 0) then 
                ! Check
                if (any(fc > 0)) then
                    ! Check if there are any internal faces to begin with
                    if (count(.not. face%BF) == 0) then 
                        print *, 'GetClosedPolygonFromTopomeshFaces: case without any ' // & 
                            'internal faces detected'

                        ! Find the next face (any next face)
                        tf = findloc( (fc > 0)  , .true., 1)
                    else
                        ! No more faces found - probably couldn't make
                        ! a closed polygon. Exit the loop
                        exit 
                    end if 
                else
                    ! All faces added, exit
                    exit 
                end if 
            end if 
            
            ! Do not subtract a counter - we need to end up in this face again.
            ! Also, don't add, we do this later on
            
            ! Set starting face for this cell
            startface = tf
                        
            ! Get neighbouring faces in correct order
            tfv = face%vert(tf, :)
            faceneig1 = GetTMVertFaceNeig(vert, tfv(1), tf)
            faceneig2 = GetTMVertFaceNeig(vert, tfv(2), tf)
            
            ! Sanity checks
            if ((size(faceneig1) == 0) .or. (size(faceneig2) == 0)) then 
                ! No neighbours found
                call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                    'could not find neighbouring faces, something wrong ' //& 
                    'with topological mesh construction. Check input')
            end if 
            if (all(faceneig1 == tf) .and. all(faceneig2 == tf)) then 
                ! Isolated face
                call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: isolated ' // & 
                    'face found, check input')
            end if
            
            ! Set the starting vertex index (if boundary vertices, needs checks)
            tvind = 0
            if ((.not. all(faceneig1 == tf)) .and. (.not. vert%BV(tfv(1)))) then 
                tvind = 1
            elseif ((.not. all(faceneig2 == tf)) .and. (.not. vert%BV(tfv(2)))) then 
                tvind = 2
            end if 

            ! If none found,  check if the first vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig1) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig1(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig1(1) == faceneig1(2)) then 
                    tvind = 1
                else
                    ! Check if any face neighbour can be taken 
                    if (any((fc(faceneig1) > 0))) then 
                        tvind = 1
                    end if 
                end if 
            end if 

            ! If none found, check if the second vertex has neighbours 
            ! with available faces
            if (any(fc(faceneig2) > 0) .and. (tvind == 0)) then 
                ! Check if both neighours are the same - in that case we can
                ! safely take this vertex as next vertex
                if (faceneig2(1) == tf) then 
                    ! Vertex has single face here, do not take it
                elseif (faceneig2(1) == faceneig2(2)) then 
                    tvind = 2
                else
                    ! Check if any face neighbour can be taken
                    if (any((fc(faceneig2) > 0))) then  
                        tvind = 2
                    end if 
                end if 
            end if 

            ! If none found, check if one of the vertices has a single non-neighbour
            ! face to continue
            if (tvind == 0) then 
                vertface = vert%GetFace(tfv(1))
                vertface = pack(vertface, (vertface /= tf .and. considerface(vertface)))
                if (size(vertface) == 1) then 
                    tvind = 1
                end if 
            end if 
            if (tvind == 0) then 
                vertface = vert%GetFace(tfv(2))
                vertface = pack(vertface, (vertface /= tf .and. considerface(vertface)))
                if (size(vertface) == 1) then 
                    tvind = 2
                end if 
            end if 
            
            ! Sanity check
            if (tvind == 0) then 
                ! No starting vertex found - no polygon can be formed since
                ! the face doesn't seem to be adjacent to other given faces
                exit
            end if
            
            ! Set current vertex
            tv = face%vert(tf, tvind)
            startvert = tv
            starttvind = tvind
        
            ! Loop 
            !=====
            do while (.true.)

                ! First, check if we can simply continue without needing
                ! to check the turndirection
                vertface = vert%GetFace(tv)
                vertface = pack(vertface, (vertface /= tf .and. &
                    considerface(vertface)) .and. (fc(vertface) > 0))
                if (size(vertface) == 1) then 
                    ! Face found, add and skip the rest of the loop
                    nf = vertface(1)
                    nfv = face%vert(nf, :)
                    istfv = nfv == tv
                    if (.not.any(istfv)) then 
                        ! Current vertex is not found in the next face, this should
                        ! not be possible
                        call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                            'next face does not have current vertex, ' // & 
                            'check input')
                    end if 
                    if (all(istfv)) then 
                        ! Next face is a face that starts and ends in the same
                        ! vertex - not supported
                        call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                            'face detected with same start and end vertex, ' // & 
                            'not supported')
                    end if 
                    
                    ! Add and update
                    tcf = [tcf, nf]
                    tcv = [tcv, tv]
                    tf = nf
                    if (.not. istfv(1)) then 
                        tv = nfv(1)
                    else
                        tv = nfv(2)
                    end if 
                    
                    ! Update counter
                    fc(tf) = fc(tf) - 1

                    ! Check exit condition
                    if ((startface == nf) .and. (tv == startvert)) then 
                        exit
                    end if 

                    ! Skip rest
                    cycle
                end if 
                
                ! Get neighbouring faces in correct order
                faceneig = GetTMVertFaceNeig(vert, tv, tf)
                        
                ! Sanity checks
                if (size(faceneig) == 0) then 
                    ! No neighbours found
                    call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: could ' // & 
                        'not find neighbouring faces, something wrong ' // & 
                        'with topological mesh construction. Check input')
                end if 
                if (all(fc(faceneig) <= 0)) then 
                    ! No neighbours with counter left
                    call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: all ' // & 
                        'neighbouring faces cannot be taken anymore, ' // &
                        'faces do not seem to form cell')
                end if 
                
                ! Find the next face
                ! We have a turn direction, so we can only check if we should
                ! throw errors
                nf = faceneig(turndirection)
                
                ! Check counter
                if (fc(nf) <= 0) then 
                    ! Cannot continue here, exit
                    exit
                end if 
                
                ! If we passed this, we should check the vertices
                nfv = face%vert(nf, :)
                istfv = nfv == tv
                if (.not.any(istfv)) then 
                    ! Current vertex is not found in the next face, this should
                    ! not be possible
                    call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                        'next face does not have current vertex, ' // & 
                        'check input')
                end if 
                if (all(istfv)) then 
                    ! Next face is a face that starts and ends in the same
                    ! vertex - not supported
                    call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                        'face detected with same start and end vertex, ' // & 
                        'not supported')
                end if 
                
                ! Add and update
                tcf = [tcf, nf]
                tcv = [tcv, tv]
                tf = nf
                if (.not. istfv(1)) then 
                    tv = nfv(1)
                else
                    tv = nfv(2)
                end if 
                
                ! Update counter
                fc(tf) = fc(tf) - 1
                
                ! Set that this direction can't be turned in anymore
                ! from neither side for this face
                if ((startface == nf) .and. (tv == startvert)) then 
                    ! Just break, turning direction etc already adjusted
                    ! before.
                    exit
                else
                    if (turndirection == 1) then 
                        if (face%vert(tf, 1) == tv) then 
                            hasturned1(tf, 1) = .true. 
                            hasturned2(tf, 2) = .true. 
                        else 
                            hasturned1(tf, 2) = .true. 
                            hasturned2(tf, 1) = .true. 
                        end if 
                    elseif (turndirection == 2) then 
                        if (face%vert(tf, 1) == tv) then 
                            hasturned1(tf, 2) = .true. 
                            hasturned2(tf, 1) = .true. 
                        else 
                            hasturned1(tf, 1) = .true. 
                            hasturned2(tf, 2) = .true. 
                        end if 
                    else
                        call gdErrorHandler('GetClosedPolygonFromTopomeshFaces: ' // & 
                            'bug detecetd when adjusting turning direction')
                    end if 
                end if 
            end do 
            
            ! Add found cell to the structure
            cc = cc + 1
            thiscellvert = ConstructIntegerDynamicArray(tcv)
            thiscellface = ConstructIntegerDynamicArray(tcf)
            cellvert = [cellvert, thiscellvert]
            cellface = [cellface, thiscellface]

            ! Housekeeping
            deallocate(tcv, tcf)
        end do 

        ! Check output
        !=============
        ! If nothing found, just return
        if (size(cellface) == 0) then 
            return
        end if 

        ! Check if we find an entry with all faces included
        do i = 1, size(cellface)
            ! Get current cells
            tempf = cellface(i)%Get()

            ! Check if they are all in faceID
            call SetDiff(faceID, tempf, diffc)
            if (size(diffc) == 0) then 
                sortedfaces = tempf 
                exit
            end if
        end do  

        ! Housekeeping
        !=============
        end associate

    end function

    !------------------------------------------------------------------!
    !                            AUXILIARY                             !
    !------------------------------------------------------------------!
    ! Tiny newton solver
    subroutine TinyNewtonSolver(x, y, converged, x0, y0, magneticField)

        ! Description
        !============
        ! Tiny newton solver implementation to refine extrema. 
        ! No guarantee on convergence (also not strictly needed)

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(out)               :: x, y
        real(R8), intent(in)                :: x0, y0 
        logical, intent(out)                :: converged 
        type(magneticFieldUDT)              :: magneticField    
        
        ! Auxiliary
        integer(I8), parameter              :: maxit = 10
        real(R8), parameter                 :: tol = 1e-10 
        integer(I8)                         :: flag
        real(R8)                            :: xi(1), yi(1), dFdxi(1), &
            dFdyi(1), d2Fdx2(1), d2Fdy2(1), d2Fdxdy(1), A(1:2, 1:2), &
            res(1:2), d(1:2)

        ! Loop 
        integer(I8)                         :: it 

        ! Initialize
        !===========
        converged = .false.
        it = 1
        xi = x0
        yi = y0
        do while ((it <= maxit) .and. (.not. converged) )
            ! Compute residuals
            call magneticField%interp%Evaluate(xi, yi, 1, 0, dFdxi)
            call magneticField%interp%Evaluate(xi, yi, 1, 0, dFdyi)

            ! Check convergence 
            res = -[dFdxi, dFdyi]
            if (all(abs(res) <= tol)) then 
                converged = .true.
                exit 
            end if 
            if (any(isnan(res))) then 
                ! Probably out of bounds of the interpolator, exit
                converged = .false. 
                exit 
            end if 
            
            ! Compute update
            call magneticField%interp%Evaluate(xi, yi, 2, 0, d2Fdx2)
            call magneticField%interp%Evaluate(xi, yi, 1, 1, d2Fdxdy)
            call magneticField%interp%Evaluate(xi, yi, 0, 2, d2Fdy2)
            A(:, 1) = [d2Fdx2, d2Fdxdy] 
            A(:, 2) = [d2Fdxdy, d2Fdy2]
            call SolveDenseLinearSystemDI(A, res, d, flag)

            ! Check if solver converged, otherwise exit 
            if (flag /= 0) then 
                exit 
            end if 
            
            ! Apply updates
            xi = xi + d(1);
            yi = yi + d(2);
            
            ! Update counter
            it = it + 1
        end do 

        ! Unpack solution
        x = xi(1)
        y = yi(1)

    end subroutine 

    ! Eigenvalue computation for 2-by-2 matrices 
    function ComputeEigenvaluesSymmetric2by2Matrix(a, b, c) result(eig)
        
        ! Description
        !============
        ! Compute the eigenvalues of a 2-by-2 symmetric matrix, of which 
        ! the diagonal values should be given as the two first entries 
        ! and the off-diagonal value as the third entry. The 
        ! eigenvalue should be per definition real, since the matrix
        ! is symmetric. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: a, b, c
        real(R8)                    :: eig(1:2)

        ! Auxiliary
        real(R8)                    :: det 

        ! Compute
        !========
        ! Determinant
        det = sqrt((a+b)**2 - 4*(a*b-c**2))

        ! First eigenvalue
        eig(1) = ((a+b) + det)/2.0_R8
        eig(2) = ((a+b) - det)/2.0_R8

    end function 

    ! Coordinate trimmer
    subroutine DeleteCurveSegment(x, y, sr, loc, dloffset, keepstartpoint, &
        keependpoint)

        ! Description
        !============
        ! An auxiliary routine that deletes a segment of a curve, defined
        ! by its x, y coordinates. The segment that is deleted is defined 
        ! by 'loc', which can be 'start', 'end', or 'both'. 'sr' is then the
        ! length coordinate (curve length goes from 0 to number of segments) at which the 
        ! cutoff is defined. Typically, this length originates from an 
        ! intersection routine. an 'offset' can be defined (0 <= offset <= 1) 
        ! which gives the relative distance that should be kept between the
        ! intersection position and the next vertex. If 'loc' is 'both', 
        ! then offset and sr is expected to be a 1-by-2
        ! array. Finally, the boolean 'keependpoints' determines whether
        ! the start and end of the curve should be conserved or deleted.

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, dimension(:), intent(inout)  :: x, y
        real(R8), dimension(:), intent(in)                  :: sr 
        character(*), intent(in)                            :: loc 
        real(R8), intent(in)                                :: dloffset(:) 
        logical, intent(in)                                 :: keepstartpoint 
        logical, intent(in)                                 :: keependpoint 

        ! Auxiliary
        real(R8), dimension(:), allocatable     :: offset
        real(R8), allocatable, dimension(:)     :: actualsr, &
            dlc, xint, yint
        integer(I8)                             :: ns
        integer(I8), allocatable, dimension(:)  :: keepind

        ! Loop
        integer(I8)                             :: i

        ! Delete segment
        !===============
        ! Initialize
        actualsr = sr 

        ! Project offset within bounds if necessary
        offset = dloffset
        !where (dloffset >= 1.0_R8) offset = 1.0_R8
        !where (dloffset <= 0.0_R8) offset = 0.0_R8
        !where (dloffset < 1.0_R8 .and. dloffset > 0.0_R8) offset = dloffset

        ! Compute length distribution
        dlc = [(i, i = 0, size(x)-1)]
        ns = size(dlc)-1

        ! Check case
        select case (loc)

        case ('start')

            ! Check sr
            if (sr(1) == real(ceiling(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) + offset(1)
            else
                if (offset(1) == 1.0_R8) then
                    actualsr(1) = real(ceiling(sr(1)), kind=R8)
                else
                    actualsr(1) = sr(1) + (ceiling(sr(1)) - sr(1))*offset(1)
                end if 
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr(1:1), xint, dlc, x)
            call Interpolate1D(actualsr(1:1), yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = ceiling(actualsr(1))+1, ns+1)]

            ! Compute new coordinates
            if (keepstartpoint) then 
                x = [x(1), xint(1), x(keepind)]
                y = [y(1), yint(1), y(keepind)]
            else
                x = [xint(1), x(keepind)]
                y = [yint(1), y(keepind)]
            end if 

        case ('end')

            ! Check sr
            if (sr(1) == real(floor(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) - offset(1)
            else
                actualsr(1) = sr(1) - (sr(1) - floor(sr(1)))*offset(1)
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr(1:1), xint, dlc, x)
            call Interpolate1D(actualsr(1:1), yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = 1, floor(actualsr(1))+1)]

            ! Compute new coordinates
            if (keependpoint) then 
                x = [x(keepind), xint, x(ns+1)]
                y = [y(keepind), yint, y(ns+1)]
            else
                x = [x(keepind), xint]
                y = [y(keepind), yint]
            end if 

        case ('both')

            ! Check sr
            if (sr(1) == real(ceiling(sr(1)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(1) = sr(1) + offset(1)
            else
                if (offset(1) == 1.0_R8) then
                    actualsr(1) = real(ceiling(sr(1)), kind=R8)
                else
                    actualsr(1) = sr(1) + (ceiling(sr(1)) - sr(1))*offset(1)
                end if 
            end if
            if (sr(2) == real(floor(sr(2)), kind=R8)) then 
                ! Intersection exactly at end node, set segment index
                actualsr(2) = sr(2) - offset(2)
            else
                actualsr(2) = sr(2) - (sr(2) - floor(sr(2)))*offset(2)
            end if 

            ! Compute intersection coordinates
            call Interpolate1D(actualsr, xint, dlc, x)
            call Interpolate1D(actualsr, yint, dlc, y)

            ! Compute which points to keep 
            keepind = [(i, i = ceiling(actualsr(1))+1, floor(actualsr(2))+1)]

            ! Compute new coordinates
            if (keepstartpoint .and. keependpoint) then 
                x = [x(1), xint(1), x(keepind), xint(2), x(ns+1)]
                y = [y(1), yint(1), y(keepind), yint(2), y(ns+1)]
            elseif (keepstartpoint) then 
                x = [x(1), xint(1), x(keepind), xint(2)]
                y = [y(1), yint(1), y(keepind), yint(2)]
            elseif (keependpoint) then 
                x = [xint(1), x(keepind), xint(2), x(ns+1)]
                y = [yint(1), y(keepind), yint(2), y(ns+1)]
            else
                x = [xint(1), x(keepind), xint(2)]
                y = [yint(1), y(keepind), yint(2)]
            end if 

        case default 

            call gdErrorHandler('DetermineCurveSegment: unknown location')

        end select

    end subroutine

    ! Tracer updater
    subroutine UpdateTracersFromTopomesh(topomesh, tracer, &
        magneticField, vessel, options)

        ! Description
        !============
        ! This routine updates the field/vessel tracers by adding the 
        ! topological mesh points as saddle points 

        ! Modules
        !========
        use mod_structured2Dgridding

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshUDT), intent(in)          :: topomesh
        class(ContourTracerUDT), allocatable, intent(inout)  :: tracer
        type(TopomeshOptionsUDT), intent(in)    :: options
        type(MagneticFieldUDT), intent(in)      :: magneticField
        type(VesselUDT), intent(in)             :: vessel

        ! Auxiliary
        real(R8)                                :: dxfracmin, dyfracmin
        real(R8), allocatable, dimension(:)     :: xb, yb, &
            xg, yg, Vf, xgv, ygv, xps, yps
        real(R8), parameter                     :: emptyR8(0)= 0
        real(R8), allocatable, dimension(:)     :: xtp, ytp, Ftp
        integer(I8)                             :: ntp 
        integer(I8), allocatable, dimension(:)  :: IDs
        integer(I8), parameter                  :: emptyI8(0) = 0
        logical, allocatable, dimension(:)      :: includevert

        ! Grid
        !=====
        ! Initialize
        dxfracmin = 1e-4_R8 
        dyfracmin = 1e-4_R8

        ! Determine domain bounds based on vessel and magnetic field extent
        call vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv]), maxval([xps, magneticField%interp%xgv])]
        yb = [minval([yps, magneticField%interp%ygv]), maxval([yps, magneticField%interp%ygv])]

        ! Construct refined grid based on tangency points and extrema
        includevert = (topomesh%vert%type == TMvertextp1ID) .or. &
            (topomesh%vert%type == TMvertextp2ID) .or. (topomesh%vert%type == TMvertexsaddleID)
        ! includevert = topomesh%vert%type == TMvertexsaddleID
        ntp = count(includevert)
        allocate(xtp(ntp), ytp(ntp), Ftp(ntp), IDs(ntp))
        xtp = pack(topomesh%vert%x, includevert)
        ytp = pack(topomesh%vert%y, includevert)
        Ftp = pack(topomesh%vert%fval, includevert)
        IDs = pack(topomesh%vert%ID, includevert)
        call ConstructRefined2DStructuredGrid(xg, yg, xgv, ygv, xb, yb, &
            options%vresx, options%vresy, xtp, ytp, 5, 5, dxfracmin, dyfracmin)  

        ! Evaluate magnetic field and vessel
        allocate(Vf(size(xg)))
        call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
        call magneticField%interp%Evaluate(xtp, ytp, 0, 0, Ftp)

        ! Update the tracer 
        tracer = ConstructStructuredTracer(&
            reshape(Vf, [size(xgv), size(ygv)]), xgv, ygv, &
            xtp, ytp, Ftp, IDs, topomesh%vert%type(IDs) ==TMvertexsaddleID, &
            tracer%npmin, tracer%npmax, tracer%dl)

    end subroutine
    !------------------------------------------------------------------!
    !                               I/O                                !
    !------------------------------------------------------------------!
    ! Topological mesh data writer
    subroutine WriteTopologicalMesh(topomesh, filename, isuptodateopt)

        ! Description
        !============
        ! Routine to write topological mesh file in our own format. Also
        ! handy for reading in once constructed. The following data is
        ! written in a column-wise fashion:

        ! Vertex data
        !============
        ! 'vertices' 
        ! <vert%ntot> 
        ! 'ID, x, y, type, fsID, fval, BV'
        ! <ID, x, y, type, fsID, fval, BV (as zero or one)>
        ! 'vertex face pointer'
        ! <vertex%faceP(:, 1), vertex%faceP(:, 2)>
        ! 'vertex faces'
        ! <vertex%face>
        ! 'faces'
        ! <face%ntot> 
        ! 'ID, fsID, type, vert1, vert2, BF, nc'
        ! <ID, fsID, type, vert(:, 1), vert(:, 2), BF, face%x%size(), label>
        ! 'face <nf>' <repeated for each face including header>
        ! <x, y> 
        ! 'cells'
        ! <cell%ntot,  cells%nvert, cells%nface> 
        ! 'ID, cell polygon size>'
        ! <ID, nc>
        ! 'cell vertices'
        ! <cell%vert>
        ! 'cell vertex pointer'
        ! <cell%vertP>
        ! 'cell faces'
        ! <cell%face> 
        ! 'cell face pointer'
        ! <cell%faceP(:, 1), cell%faceP(:, 2)>
        ! 'cell polygon <nc>' (repeated nc times)
        ! <x, y> 
        ! 'tubes'
        ! <tubes%ntot,  tubes%nface, tubes%ncell> 
        ! 'tube faces'
        ! <tube%face>
        ! 'tube face pointer'
        ! <tube%faceP(:, 1), tube%faceP(:, 2)>
        ! 'tube cells'
        ! <tube%cell>
        ! 'tube cell pointer'
        ! <tube%cellP(:, 1), tube%cellP(:, 2)>
        
        ! Declare variables
        !==================
        ! Modules 
        use mod_plotter 
        use mod_specialchars, only  : filesepchar
        use mod_definitions, only   : goatversion

        ! Arguments
        class(TopomeshUDT)                      :: topomesh
        character(*), intent(in)                :: filename 
        logical, intent(in), optional           :: isuptodateopt

        ! Auxiliary
        integer                                 :: fu, BVval, cvsize, &
            cfsize
        integer(I8), allocatable, dimension(:)  :: tID
        real(R8), allocatable, dimension(:)     :: xf, yf, tfval
        character(:), allocatable               :: dir
        logical                                 :: isuptodate
        logical, allocatable, dimension(:)      :: BV, BF
        type(RealDynamicArrayUDT), allocatable  :: xcda(:), ycda(:)

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell,   &
            t       => topomesh%tube    &
        )

        ! Switch to check if full topomesh is up to date. If not, 
        ! certain (many) things will not be written out (cells, 
        ! some connectivity data, ...)
        if (present(isuptodateopt)) then 
            isuptodate = isuptodateopt
        else
            isuptodate = .true.
        end if 

        ! Construct writing directory
        dir = plotdir // filesepchar // filename // '.dat'

        ! Open file
        open (action='write', file=trim(dir), newunit=fu, &
             status='unknown')

        ! Write header
        write(fu, *) 'VERSION' // goatversion

        ! Check data
        !===========
        ! Only the fields that are not added when constructing vertices,
        ! faces or cells. The rest should be always available
        if (.not. allocated(v%BV)) then 
            allocate(BV(v%ntot))
            BV = .false.
        elseif (size(v%BV) /= v%ntot) then 
            allocate(BV(v%ntot))
            BV = .false.
        else
            BV = v%BV
        end if 
        if (.not. allocated(f%BF)) then 
            allocate(BF(f%ntot))
            BF = .false.
        elseif (size(f%BF) /= f%ntot) then 
            allocate(BF(f%ntot))
            BF = .false.
        else 
            BF = f%BF
        end if

        ! Write vertex data
        !==================
        ! Number of vertices
        write (fu, *) 'vertices'
        write (fu, *) v%ntot

        ! Basic vertex data
        write (fu, *) 'ID, x, y, type, fsID, fval, BV'
        do i = 1, v%ntot
            if (BV(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if
            write (fu, *) v%ID(i), v%x(i), v%y(i), v%type(i), v%fsID(i), v%fval(i), BVval
        end do 

        ! Vertex face data
        if (allocated(v%faceP) .and. allocated(v%face)) then 
            if (size(v%faceP, 1) == v%ntot) then ! may not be up to date
                write (fu, *) 'vertex face pointer'
                do i = 1, v%ntot
                    write (fu, *) v%faceP(i, 1), v%faceP(i, 2)
                end do 

                write (fu, *) 'vertex facelist'
                do i = 1, size(v%face)
                    write (fu, *) v%face(i)
                end do 
            else
                isuptodate = .false.
            end if 
        else
            isuptodate = .false.
        end if 

        ! Write face data
        !================
        ! Number of faces
        write (fu, *) 'faces'
        write (fu, *) f%ntot 

        ! Basic face data
        write (fu, *) 'ID, fsID, type, vert1, vert2, BF, nc'
        do i = 1, f%ntot 
            if (BF(i)) then 
                BVval = 1
            else 
                BVval = 0
            end if 
            write (fu, *) f%ID(i), f%fsID(i), f%type(i), f%vert(i, 1), &
                f%vert(i, 2), BVval, f%x(i)%size()
        end do 
        
        ! Face coordinates
        write (fu, *) 'face coordinates'
        do i = 1, f%ntot 
            ! Write header
            write (fu, *) 'face ', i 

            ! Write coordinates
            xf = f%x(i)%Get()
            yf = f%y(i)%Get()
            do j = 1, size(xf)
                write(fu, *) xf(j), yf(j)
            end do 
        end do 

        ! Write cell data
        !================
        ! Number of cells, number of cell vertices, number of cell faces
        if (allocated(c%vert) .and. isuptodate) then 
            if ((size(c%vert) == c%ntotv) .and. .not. (any(c%vert == 0))) then 
                cvsize = size(c%vert)
            else
                cvsize = 0
                isuptodate = .false.
            end if 
        else
            cvsize = 0
            isuptodate = .false.
        end if 
        if (allocated(c%face) .and. isuptodate) then 
            if (size(c%face) == c%ntotf) then 
                cfsize = size(c%face)
            else
                cfsize = 0
                isuptodate = .false.
            end if 
        else
            cfsize = 0
            isuptodate = .false.
        end if 
        write (fu, *) 'cells'
        if (isuptodate) then 
            write (fu, *) c%ntot, cvsize, cfsize
        else
            write (fu, *) 0_I8, cvsize, cfsize
        end if 

        ! Print additional data
        if (isuptodate) then 
            ! Compute cell polygon data
            allocate(xcda(c%ntot), ycda(c%ntot))
            do i = 1, c%ntot
                call ConstructTopologicalMeshCellPolygon(topomesh, i, xcda(i), ycda(i))
            end do 

            ! Cell data
            write (fu, *) 'ID, nc'
            do i = 1, c%ntot 
                write (fu, *) i, xcda(i)%Size()
            end do 

            ! Cell vertices
            write (fu, *) 'cell vertexlist'
            do i = 1, cvsize
                write (fu, *) c%vert(i)
            end do 
            write (fu, *) 'cell vertex pointer'
            do i = 1, c%ntot 
                write (fu, *) c%vertP(i, 1), c%vertP(i, 2)
            end do 

            ! Cell faces
            write (fu, *) 'cell facelist'
            do i = 1, cfsize
                write (fu, *) c%face(i)
            end do 
            write (fu, *) 'cell face pointer'
            do i = 1, c%ntot 
                write (fu, *) c%faceP(i, 1), c%faceP(i, 2)
            end do 

            ! Cell polygons
            write (fu, *) 'cell polygons'
            do i = 1, c%ntot 
                ! Write header
                write (fu, *) 'cell ', i 

                ! Write coordinates
                xf = xcda(i)%Get()
                yf = ycda(i)%Get()
                do j = 1, size(xf)
                    write(fu, *) xf(j), yf(j)
                end do 
            end do 
        end if 

        ! Write flux surface data
        !========================
        if (isuptodate) then 
            ! Header and sizes
            write (fu, *) 'flux surfaces'
            write (fu, *) topomesh%nfs 
            
            ! Data
            write (fu, *) 'ID, fval'
            tID = topomesh%fsID%Get()
            tfval = topomesh%fsfval%Get()
            do i = 1, size(tID)
                write (fu, *) tID(i), tfval(i)
            end do 
        else
            ! Header and sizes
            write (fu, *) 'flux surfaces'
            write (fu, *) 0_I8 
        end if 

        ! Write tube data
        !================
        if (isuptodate) then 
            ! Header and sizes
            write (fu, *) 'tubes'
            write (fu, *) t%ntot, t%nface, t%ncell 

            ! Tube faces
            write (fu, *) 'tube facelist'
            do i = 1, t%nface
                write (fu, *) t%face(i)
            end do 
            write (fu, *) 'tube face pointer'
            do i = 1, t%ntot 
                write (fu, *) t%faceP(i, 1), t%faceP(i, 2)
            end do 

            ! Tube cells
            write (fu, *) 'tube cells'
            do i = 1, t%ncell
                write (fu, *) t%cell(i)
            end do 
            write (fu, *) 'tube cell pointer'
            do i = 1, t%ntot 
                write (fu, *) t%cellP(i, 1), t%cellP(i, 2)
            end do 
        else
            ! Header and sizes
            write (fu, *) 'tubes'
            write (fu, *) 0_I8, 0_I8, 0_I8
        end if 

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)

    end subroutine

    ! Topological mesh data reader 
    subroutine ReadTopologicalMesh(topomesh, filepath)

        ! Description
        !============
        ! This routine reads in the topological mesh data of an 
        ! existing topological mesh (which is assumed to be correctly
        ! set up) and adds any additional interconnection data that is 
        ! not directly available in the file itself.

        ! Declare variables
        !==================
        ! Modules 
        use mod_readwrite
        use mod_inputfileparser

        ! Arguments
        class(TopomeshUDT), intent(out)         :: topomesh
        character(*), intent(in)                :: filepath 

        ! Auxiliary
        integer                                 :: fu, BVval, &
            nfc, ncc
        real(R8), allocatable, dimension(:)     :: xf, yf, tfval 
        integer(I8), allocatable, dimension(:)  :: tID
        character(:), allocatable               :: thisline
        logical                                 :: reachedeof

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Unpack
        associate(&
            v       => topomesh%vert,   &
            f       => topomesh%face,   &
            c       => topomesh%cell,   &
            t       => topomesh%tube    &
        )

        ! Open file
        open (action='read', file=trim(filepath), newunit=fu, &
             status='old')

        ! Initialize
        call topomesh%Initialize()

        ! Vertices
        !=========
        ! Read until header found
        call ReadUntilFound(fu, 'vertices', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find vertices')
        end if 

        ! Read number of vertices
        read(fu, *) v%ntot

        ! Initialize vertices 
        call v%Initialize(nv=v%ntot)

        ! Basic vertex data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, v%ntot
            read(fu, *) v%ID(i), v%x(i), v%y(i), v%type(i), v%fsID(i), v%fval(i), BVval
            if (BVval == 1) then 
                v%BV(i) = .true.
            else
                v%BV(i) = .false.
            end if 
        end do 

        ! Vertex face data
        call ReadSingleLine(fu, thisline, reachedeof)
        if (allocated(v%faceP)) then 
            deallocate(v%faceP)
        end if
        allocate(v%faceP(v%ntot, 2))
        v%faceP = 0
        do i = 1, v%ntot 
            read(fu, *) v%faceP(i, 1), v%faceP(i, 2)
        end do 
        if (allocated(v%face)) then 
            deallocate(v%face)
        end if
        allocate(v%face(sum(v%faceP(:, 2))))
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, size(v%face)
            read(fu, *) v%face(i)
        end do 

        ! Faces
        !======
        ! Read until header found
        call ReadUntilFound(fu, 'faces', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find faces')
        end if 

        ! Read number of faces
        read(fu, *) f%ntot

        ! Initialize vertices 
        call f%Initialize(nf=f%ntot)

        ! Read basic data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, f%ntot
            ! Read
            read(fu, *) f%ID(i), f%fsID(i), f%type(i), f%vert(i, 1), &
                f%vert(i, 2), BVval, nfc
            if (BVval == 1) then 
                f%BF(i) = .true.
            else
                f%BF(i) = .false.
            end if 

            ! Initialize
            f%x(i) = ConstructRealDynamicArray(spread(0.0_R8, 1, nfc))
            f%y(i) = ConstructRealDynamicArray(spread(0.0_R8, 1, nfc))
        end do 

        ! Read face coordinates
        call ReadSingleLine(fu, thisline, reachedeof) ! skip the header
        do i = 1, f%ntot 
            ! Read header
            call ReadSingleLine(fu, thisline, reachedeof) 

            ! Read coordinates
            xf = f%x(i)%Get()
            yf = f%y(i)%Get()
            do j = 1, f%x(i)%Size()
                read(fu, *) xf(j), yf(j)
            end do 
            call f%x(i)%Set(xf)
            call f%y(i)%Set(yf)
            call f%pol(i)%Construct(xf, yf)
        end do 

        ! Cells
        !======
        ! Read until header found
        call ReadUntilFound(fu, 'cells', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find cells')
        end if 

        ! Read number of faces
        read(fu, *) c%ntot, c%ntotv, c%ntotf 

        ! Initialize vertices 
        call c%Initialize(nc=c%ntot, ncv=c%ntotv, ncf=c%ntotf)

        ! Basic data
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot 
            read(fu, *) c%ID(i), ncc
        end do 

        ! Cell vertices
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntotv
            read (fu, *) c%vert(i)
        end do 

        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot
            read(fu, *) c%vertP(i, 1), c%vertP(i, 2)
        end do 

        ! Cell faces
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntotf
            read (fu, *) c%face(i)
        end do 

        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, c%ntot
            read(fu, *) c%faceP(i, 1), c%faceP(i, 2)
        end do 

        ! Flux surface data
        !==================
        ! Read until header found
        call ReadUntilFound(fu, 'flux surfaces', reachedeof)
        if (reachedeof) then 
            call gdErrorHandler('ReadTopologicalMesh: could not find flux surfaces')
        end if 

        ! Read number of surfaces
        read(fu, *) topomesh%nfs 

        ! Allocate
        allocate(tfval(topomesh%nfs), tID(topomesh%nfs))

        ! Read 
        call ReadSingleLine(fu, thisline, reachedeof)
        do i = 1, topomesh%nfs
            read(fu, *) tID(i), tfval(i)
        end do 
        topomesh%fsID = ConstructIntegerDynamicArray(tID)
        topomesh%fsfval = ConstructRealDynamicArray(tfval)

        ! Add interconnection data
        !=========================
        call AddTopologicalMeshData(topomesh)
        call AddTopologicalMeshInterconnectionData(topomesh)

        ! Housekeeping
        !=============
        ! Deallocate again

        ! Others
        end associate
        close(fu)


    end subroutine


end module
