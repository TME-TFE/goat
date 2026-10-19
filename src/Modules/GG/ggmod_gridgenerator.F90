!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides a grid generator object that allows for more
! flexible (re)gridding. Currently only 2D gridding is available, but 
! if ever need be, 3D grid generators can also be defined here. Note
! that the main grid generation implementation is defined in 
! ggmod_gridgeneration2D and ggmod_topology2D. The grid generator 
! object may contain any required state for (re)gridding or updating.
! Also, this may make adaptive refinement during simulation or 
! optimization easier (if adaptive refinement will become available
! at some point)

module ggmod_gridgenerator

    ! Modules
    !========
    use goatmod_types 
    use goatmod_userinput
    use ggmod_gridgeneration2D
    use ggmod_topology2D
    use mod_contour2D
    use mod_streamlinetracing2D
    use mod_errorhandler
    implicit none

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Grid generator object
    !======================
    ! 2D grid generator
    type :: GoatGridGenerator2DUDT

        ! Description
        !============
        ! We need to keep track of the grid generation options, the
        ! different tracer objects, ... On top of that, we track the
        ! state of the grid generator:
        ! - hastopomesh: was a topomesh constructed?
        ! - hasgrid: was a grid constructed?
        ! These states are necessary to ensure proper and efficient
        ! working of the grid generator. Proper initialization is 
        ! done in the constructor (not part of this object)

        type(GGoptionsUDT)          :: ggoptions 
        type(TopomeshOptionsUDT)    :: topomeshoptions
        type(TopomeshUDT)           :: topomesh 
        type(GridUDT)               :: grid  
        class(ContourTracerUDT), allocatable    :: fieldtracer, &
            boundarytracer 
        class(StreamlineTracerUDT), allocatable :: streamlinetracer
        type(EnvironmentUDT)        :: environment 
        type(MagneticFieldUDT)      :: magneticField 

        ! State flags
        logical :: isinitialized    = .false. 
        logical :: hasenvironment   = .false. 
        logical :: hasmagneticfield = .false. 
        logical :: hastopomesh      = .false. 
        logical :: hasgrid          = .false. 

    contains
    
        ! Constructor
        procedure :: Initialize         => InitializeGGG2D

        ! Topological mesh generation 
        procedure :: ConstructTopomesh  => ConstructTopomeshGGG2D

        ! Grid generation 
        procedure :: ConstructGrid      => ConstructGridGGG2D

        ! Update environment and boundary tracer
        procedure :: UpdateEnvironment  => UpdateEnvironmentGGG2D 

        ! Update magnetic field and tracers
        procedure :: UpdateMagneticField    => UpdateMagneticFieldGGG2D

    end type

contains 

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                     GOAT GRID GENERATOR 2D                       !
    !------------------------------------------------------------------!

    ! Constructor
    subroutine InitializeGGG2D(gridgenerator, magneticField, environment, &
        inputfilepath) 

        ! Description
        !============
        ! Main constructor for the 2D goat grid generator. Inputs 
        ! should be the options structure with already read in 
        ! options, the properly initialized magnetic field, and the 
        ! initialized environment (use eg ExtractGGData upstream)

        ! Modules
        !========
        use mod_structured2Dgridding, only: Construct2DStructuredUniformGrid

        ! Declare variables
        !==================
        ! Arguments
        class(GoatGridGenerator2DUDT)   :: gridgenerator
        type(MagneticFieldUDT)      :: magneticField 
        type(EnvironmentUDT)        :: environment
        character(*), intent(in)    :: inputfilepath

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vf, Vv, xgv, ygv, Vfx, Vfy
        real(R8), parameter                     :: emptyR8(0)= 0
        integer(I8)                             :: nv
        integer(I8), parameter                  :: emptyI8(0) = 0
        logical, parameter                      :: emptyL(0) = .false.

        ! Loop

        ! Initialize
        !===========
        ! Check flag
        if (gridgenerator%isinitialized) then 
            print *, 'InitializeGGG2D: re-initializing grid generator'
        end if 

        ! Add
        gridgenerator%environment = environment 
        gridgenerator%magneticField = magneticField

        ! Associate for ease
        associate(&
            topomeshoptions     => gridgenerator%topomeshoptions,   &
            ggoptions           => gridgenerator%ggoptions,         &
            environment         => gridgenerator%environment,       &
            magneticField       => gridgenerator%magneticField,     &
            fieldtracer         => gridgenerator%fieldtracer,       &
            boundarytracer      => gridgenerator%boundarytracer,    &
            streamlinetracer    => gridgenerator%streamlinetracer   &
            )

        ! Set topological mesh options
        topomeshoptions%inputfilepath = inputfilepath 
        call topomeshoptions%Set()

        ! Check
        if (topomeshoptions%readexistingTM) then 
            print *, 'InitialieGGG2D: setting readexistingTM to false'
            topomeshoptions%readexistingTM = .false. 
        end if 

        ! Set grid generation options
        ggoptions%inputfilepath = inputfilepath 
        call ggoptions%Set()

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
        boundarytracer = ConstructStructuredTracer(&
            reshape(Vv, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topomeshoptions%npmin, &
            topomeshoptions%npmax, topomeshoptions%dl)

        ! Orthogonal lines
        streamlinetracer = ConstructStructuredStreamlineTracer(&
            reshape(Vfx, [topomeshoptions%vresx, topomeshoptions%vresy]), &
            reshape(Vfy, [topomeshoptions%vresx, topomeshoptions%vresy]), & 
            xgv, ygv, step=ggoptions%orthtracerstep, nsteps=ggoptions%orthtracernsteps) 

        ! Set initialized to true
        gridgenerator%isinitialized     = .true. 
        gridgenerator%hasmagneticfield  = .true. 
        gridgenerator%hasenvironment    = .true.

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Topomesh construction
    subroutine ConstructTopomeshGGG2D(gridgenerator)

        ! Description
        !============
        ! Wrapper for topomesh construction. If topomesh was already 
        ! constructed, a message is shown. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatGridGenerator2DUDT)       :: gridgenerator 

        ! Construct
        !==========
        ! Checks
        if (.not. gridgenerator%isinitialized) then 
            ! Print message
            call gdErrorHandler('ConstructTopomeshGGG2D: grid generator is ' // & 
                'not yet initialized')
        end if 
        if (gridgenerator%hastopomesh) then 
            print *, 'ConstructTopomeshGGG2D: reconstructing topological mesh'
        end if 

        ! Construct
        call ConstructTopologicalMesh(&
            gridgenerator%environment%vessel, gridgenerator%magneticField, &
            gridgenerator%topomeshoptions, gridgenerator%topomesh, &
            gridgenerator%fieldtracer, gridgenerator%boundarytracer, &
            gridgenerator%streamlinetracer)

        ! Set logicals
        gridgenerator%hastopomesh   = .true. ! Topomesh is now present
        gridgenerator%hasgrid       = .false. ! grid is now outdated 

    end subroutine

    ! Grid construction
    subroutine ConstructGridGGG2D(gridgenerator)

        ! Description
        !============
        ! This is a wrapper for generating a structured 2D aligned grid
        ! using the goat grid generator. If preliminaries are not met, 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatGridGenerator2DUDT)       :: gridgenerator 

        ! Initialize
        !===========
        ! Checks
        if (.not. gridgenerator%isinitialized) then 
            ! Call error
            call gdErrorHandler('ConstructGridGGG2D: grid generator is ' // & 
                'not yet initialized')
        end if 
        if (.not. gridgenerator%hastopomesh) then 
            ! Print message
            print *, 'ConstructGridGGG2D: constructing topological mesh'

            ! Construct topomesh
            call gridgenerator%ConstructTopomesh()
        end if 
        if (gridgenerator%hasgrid) then 
            ! Print message
            print *, 'ConstructGridGGG2D: reconstructing grid'
        end if 

        ! Generate grid
        !==============
        ! Call the unstructured aligned grid generator
        call GenerateUnstructuredAlignedGrid(gridgenerator%grid, &
            gridgenerator%topomesh, gridgenerator%magneticField, &
            gridgenerator%environment%vessel, gridgenerator%fieldtracer, &
            gridgenerator%boundarytracer, gridgenerator%streamlinetracer, &
            gridgenerator%ggoptions, restrictwallcontacttotargets=&
                gridgenerator%topomeshoptions%removewidegridregions)

        ! Set logicals
        gridgenerator%hasgrid = .true. 

    end subroutine

    ! Environment update 
    subroutine UpdateEnvironmentGGG2D(gridgenerator, environment)

        ! Description
        !============
        ! Update the environment and related data (e.g. vessel 
        ! description). Note: magnetic field may remain the same 
        ! without changes, so we don't tag that as unavailable. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatGridGenerator2DUDT)       :: gridgenerator 
        type(EnvironmentUDT), intent(in)    :: environment 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vv, xgv, ygv
        real(R8), parameter                     :: emptyR8(0)= 0
        integer(I8)                             :: nv
        integer(I8), parameter                  :: emptyI8(0) = 0
        logical, parameter                      :: emptyL(0) = .false.

        ! Initialize
        !===========
        ! Checks
        if (.not. gridgenerator%isinitialized) then 
            call gdErrorHandler('UpdateEnvironmentGGG2D: grid generator ' // & 
                'not initialized')
        end if
        if (gridgenerator%hasenvironment) then 
            print *, 'UpdateEnvironmentGGG2D: updating environment'
        end if 

        ! Update
        !=======
        ! Assign environment
        gridgenerator%environment = environment 

        ! Associate for ease
        associate(&
            topomeshoptions     => gridgenerator%topomeshoptions,   &
            ggoptions           => gridgenerator%ggoptions,         &
            environment         => gridgenerator%environment,       &
            magneticField       => gridgenerator%magneticField,     &
            fieldtracer         => gridgenerator%fieldtracer,       &
            boundarytracer      => gridgenerator%boundarytracer,    &
            streamlinetracer    => gridgenerator%streamlinetracer   &
            )

        ! Determine domain bounds based on vessel and magnetic field extent
        call environment%vessel%plfvessel%ps%GetVertices(xps, yps)
        xb = [minval([xps, magneticField%interp%xgv(2:size(magneticField%interp%xgv)-1)]), &
            maxval([xps, magneticField%interp%xgv(2:size(magneticField%interp%xgv)-1)])]
        yb = [minval([yps, magneticField%interp%ygv(2:size(magneticField%interp%ygv)-1)]), &
            maxval([yps, magneticField%interp%ygv(2:size(magneticField%interp%ygv)-1)])]

        ! Construct a 2D structured grid for tracing (may be extended
        ! in the future for different grid types)
        nv = topomeshoptions%vresx*topomeshoptions%vresy
        allocate(xg(nv), yg(nv), Vv(nv), xgv(topomeshoptions%vresx), &
            ygv(topomeshoptions%vresy))
        call Construct2DStructuredUniformGrid(xg, yg, xgv, ygv, xb, yb, &
            topomeshoptions%vresx,  topomeshoptions%vresy, 0.0_R8, 0.0_R8)

        ! Evaluate the field and vessel values
        call environment%vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)

        ! Vessel contours
        boundarytracer = ConstructStructuredTracer(&
            reshape(Vv, [topomeshoptions%vresx, topomeshoptions%vresy]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topomeshoptions%npmin, &
            topomeshoptions%npmax, topomeshoptions%dl)

        ! Set logicals
        gridgenerator%hasenvironment    = .true.
        gridgenerator%hastopomesh       = .false. ! has to be updated
        gridgenerator%hasgrid           = .false. ! has to be updated

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Magnetic field update
    subroutine UpdateMagneticFieldGGG2D(gridgenerator, magneticField)

        ! Description
        !============
        ! Update the magnetic field and related data such as the 
        ! field and streamline tracers. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatGridGenerator2DUDT)       :: gridgenerator 
        type(magneticFieldUDT), intent(in)  :: magneticField  

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: xb, yb, xps, &
            yps, xg, yg, Vv, xgv, ygv, Vf, Vfx, Vfy 
        real(R8), parameter                     :: emptyR8(0)= 0
        integer(I8)                             :: nv
        integer(I8), parameter                  :: emptyI8(0) = 0
        logical, parameter                      :: emptyL(0) = .false.

        ! Initialize
        !===========
        ! Checks
        if (.not. gridgenerator%isinitialized) then 
            call gdErrorHandler('UpdateMagneticFieldGGG2D: grid generator ' // & 
                'not initialized')
        end if
        if (gridgenerator%hasmagneticfield) then 
            print *, 'UpdateMagneticFieldGGG2D: updating magnetic field'
        end if 

        ! Update
        !=======
        ! Assign environment
        gridgenerator%magneticField = magneticField 

        ! Associate for ease
        associate(&
            topomeshoptions     => gridgenerator%topomeshoptions,   &
            ggoptions           => gridgenerator%ggoptions,         &
            environment         => gridgenerator%environment,       &
            magneticField       => gridgenerator%magneticField,     &
            fieldtracer         => gridgenerator%fieldtracer,       &
            boundarytracer      => gridgenerator%boundarytracer,    &
            streamlinetracer    => gridgenerator%streamlinetracer   &
            )

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

        ! Orthogonal lines
        streamlinetracer = ConstructStructuredStreamlineTracer(&
            reshape(Vfx, [topomeshoptions%vresx, topomeshoptions%vresy]), &
            reshape(Vfy, [topomeshoptions%vresx, topomeshoptions%vresy]), & 
            xgv, ygv, step=ggoptions%orthtracerstep, nsteps=ggoptions%orthtracernsteps) 
    
        ! Set logicals
        gridgenerator%hasmagneticfield  = .true.
        gridgenerator%hastopomesh       = .false. ! has to be updated
        gridgenerator%hasgrid           = .false. ! has to be updated

        ! Housekeeping
        !=============
        end associate

    end subroutine 


end module
