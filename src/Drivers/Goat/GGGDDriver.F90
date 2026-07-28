subroutine GGGDDriver(goatoptions)

    ! Description
    !============
    ! This driver first runs the grid generator and afterwards the 
    ! grid deformation module. This is mainly intended as a testing 
    ! routine for coupling the grid generator to the grid deformation

    ! Initialize
    !===========
    ! Modules
    use goatmod_types 
    use goatmod_userinput
    use ggmod_topology2D
    use ggmod_gridgeneration2D
    use mod_contour2D, only : ContourTracerUDT, ConstructStructuredTracer, &
        ContourUDT
    use mod_streamlinetracing2D, only: StreamlineTracerUDT, &
        ConstructStructuredStreamlineTracer
    use mod_structured2Dgridding

    ! The usual
    implicit none 

    ! Declare variables
    !==================
    ! Arguments
    type(GoatoptionsUDT)        :: goatoptions 
    type(GridUDT)               :: grid 
    type(MagneticFieldUDT)      :: magneticField 
    type(EnvironmentUDT)        :: environment

    ! Other options
    type(GGoptionsUDT)          :: ggoptions
    type(GDoptionsUDT)          :: gdoptions
    type(GridOptionsUDT)        :: gridoptions 
    type(TopomeshOptionsUDT)    :: topomeshoptions
    
    ! Auxiliary
    type(TopomeshUDT)           :: topomesh
    class(ContourTracerUDT), allocatable    :: fieldtracer, vesseltracer
    class(StreamlineTracerUDT), allocatable :: streamlinetracer
    type(VesselUDT)             :: tempvessel
    type(GGTMDataUDT)           :: ggtmdata

    real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
        yps, xg, yg, Vf, Vv, xgv, ygv, Vfx, Vfy
    real(R8), parameter                     :: emptyR8(0)= 0
    integer(I8), allocatable, dimension(:)  :: facelabelsGG, facelabelsGD
    integer(I8)                             :: nv
    integer(I8), parameter                  :: emptyI8(0) = 0
    logical, parameter                      :: emptyL(0) = .false.

    ! Loop

    ! Initialize
    !===========
    ! Read and extract data
    call ExtractGGData(magneticField, environment, goatoptions)

    ! Set topological mesh options
    topomeshoptions%inputfilepath = goatoptions%inputfilepath 
    call topomeshoptions%Set()

    ! Set grid generation options
    ggoptions%inputfilepath = goatoptions%inputfilepath 
    call ggoptions%Set()

    ! Set grid deformation options
    gdoptions%inputfilepath = goatoptions%inputfilepath
    call gdoptions%Set()

    ! Set grid options
    gridoptions%inputfilepath = goatoptions%inputfilepath 
    call gridoptions%Set()

    ! Construct (initial) tracers
    !============================
    ! Determine domain bounds based on vessel and magnetic field extent
    call environment%vessel%plfvessel%ps%GetVertices(xps, yps)
    xb = [minval([xps, magneticField%interp%xgv(2:size(magneticField%interp%xgv)-1)]), &
        maxval([xps, magneticField%interp%xgv(2:size(magneticField%interp%xgv)-1)])]
    yb = [minval([yps, magneticField%interp%ygv(2:size(magneticField%interp%ygv)-1)]), &
        maxval([yps, magneticField%interp%ygv(2:size(magneticField%interp%ygv)-1)])]

    ! Construct a 2D structured grid for tracing (may be extended
    ! in the future for different grid types)
    nv = topomeshoptions%vresx*topomeshoptions%vresy
    allocate(xg(nv), yg(nv), Vf(nv), Vv(nv), xgv(topomeshoptions%vresx), &
        ygv(topomeshoptions%vresy))
    call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
        topomeshoptions%vresx,  topomeshoptions%vresy, 0.0_R8, 0.0_R8)

    ! Evaluate the field and vessel values
    allocate(Vfx(nv), Vfy(nv))
    call environment%vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)
    call magneticField%interp%Evaluate(xg, yg, 0, 0, Vf)
    call magneticField%interp%Evaluate(xg, yg, 1, 0, Vfx)
    call magneticField%interp%Evaluate(xg, yg, 0, 1, Vfy)

    ! Field contours
    fieldtracer = ConstructStructuredTracer(&
        reshape(Vf, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
        emptyR8, emptyR8, emptyR8, emptyi8, emptyL, topomeshoptions%npmin, &
        topomeshoptions%npmax, topomeshoptions%dl)

    ! Vessel contours
    vesseltracer = ConstructStructuredTracer(&
        reshape(Vv, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
        emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topomeshoptions%npmin, &
        topomeshoptions%npmax, topomeshoptions%dl)

    ! Orthogonal lines
    streamlinetracer = ConstructStructuredStreamlineTracer(&
        reshape(Vfx, [topomeshoptions%vresx, topomeshoptions%vresy]), &
        reshape(Vfy, [topomeshoptions%vresx, topomeshoptions%vresy]), & 
        xgv, ygv, step=ggoptions%orthtracerstep, nsteps=ggoptions%orthtracernsteps) 

    ! Visualize by tracing contours
    !resc = 100
    !dv = (maxval(Vf) - minval(Vf))
    !cgv = [(k, k = 0, resc)]*(dv*0.90_R8)/real(resc, kind=R8) + minval(Vf) + dv*0.05
    !contours = fieldtracer%TraceContours(cgv)
    !allocate(pcontours(size(contours)))
    !do k = 1, size(contours)
    !    call pcontours(k)%Construct(contours(k)%x, contours(k)%y)
    !end do 
    !call tempps%Construct(pcontours)
    !call tempps%WriteData('mfcontours')


    ! Generate the topological mesh
    !==============================
    tempvessel = environment%vessel
    call ConstructTopologicalMesh(tempvessel, magneticField, &
        topomeshoptions, topomesh, fieldtracer, vesseltracer, streamlinetracer)

    ! Generate the grid
    !==================
    call GenerateUnstructuredAlignedGrid(grid, topomesh, magneticField, &
        environment%vessel, fieldtracer, vesseltracer, streamlinetracer, &
        ggoptions, ggtmdataopt=ggtmdata)

    ! Apply grid deformation
    !=======================
    ! Get label mapping
    call GetGridFaceLabelMappingGD(grid, topomesh, facelabelsGG, facelabelsGD)

    ! Overwrite grid options
    gridoptions%facelabelmappingGG = facelabelsGG 
    gridoptions%facelabelmappingGD = facelabelsGD 
    gridoptions%facelabelsubfrom   = emptyI8 
    gridoptions%facelabelsubto     = emptyI8 
    print *, 'GGGDDriver: facelabelsGG: ', facelabelsGG 
    print *, 'GGGDDriver: facelabelsGG: ', facelabelsGD

    ! Recompute topological data from grid for new face labels
    call ComputeTopologicalData(grid, topomesh, ggoptions)

    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! Extract the required grid data for grid deformation
    call ExtractGridData(grid, 'traduitb2us', gridoptions)

    ! Run deformation
    call RunGridOptimization(grid, magneticField, environment, &
        gdoptions)

    ! Write data
    !===========
    ! Translate labels etc
    call TranslateGridLabels(grid, topomesh, environment%vessel, ggtmdata, &
        ggoptions, 'solps')

    ! Recompute topological data from grid for new face labels
    call ComputeTopologicalData(grid, topomesh, ggoptions)

    ! Grid data
    call WriteGOAT(goatoptions, grid, magneticField, environment)

    ! b2ag file
    !call Writeb2agdat(goatoptions, grid)



end subroutine
