!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the constraint classes specific for the grid
! deformation. It relies on user input defined in the gdmod_userinput
! module, and on the gdmod_user types.

! The constraints are structured as follows:
! - All derived constraint types inherit from the 'mother' type
!   'GenericConstraintsGDUDT', which contains the field 'ncon', the
!   number of constraints, and the initialization and evaluation
!   routines.
! - The overarching constraints structure contains an object 'eqcon' and
!   'ineqcon', which are objects that contain the specific equality and
!   inequality constraints, respectively. Both objects have general
!   initialization, evaluation, and destruction routines that should be
!   used in the optimizer.
! - The equality constraints contain different fields (e.g. fluxfunction
!   ) with a logical (e.g. dofluxfunction) that indicates whether the
!   constraint should be considered.
! - Each specific constraints (e.g. fluxfunction) has its own evaluation
!   , initialization, and destruction routines (type-bound).

module gdmod_constraints

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_utility
    use mod_sparseinterface
    use optmod_constraints
    use gdmod_types
    use gdmod_designvariables
    use gdmod_utility_optimization
    use gdmod_plots
    use gdmod_userinput
    use PolygonLevelsetFunction2D
    use, intrinsic :: ieee_arithmetic

    ! The usual
    implicit none
    save
    public

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Auxiliary types
    !================
    ! Structures for flux surface constraints
    type :: FFCStructureUDT

        ! Description
        !============
        ! UDT to contain flux surface data for constraint evaluation.
        ! Decoupled from grid description to increase modularity.
        ! Fields:
        !   - nID   : number of vertices (or vertex IDs) present in the
        !           flux surface
        !   - ID    : vertex ID array (1D array of size nID)
        !   - PsiD  : desired psi value for this flux surface
        !   - fsID  : (optional) flux surface ID of flux surface that
        !           the structure belongs to (scalar)

        integer(I8)                     :: nID, fsID
        integer(I8), allocatable        :: ID(:)
        real(R8)                        :: PsiD

    end type

    ! Structure for keeping track of degrees of freedom
    type :: DOFGStructureUDT

        ! Description
        !============
        ! UDT that contains information (for one degree of freedom group) of
        ! how many degrees of freedom there are (in the 'dofs' field)
        ! and to which constraints (in the 'cons' field) they can be
        ! attributed. Additionally has the 'vert' field, which is only
        ! used for debugging/visualization/diagnostics.

        integer(I8)                 :: dofs
        integer(I8), allocatable    :: cons(:), vert(:)

    end type

    ! Structure to keep track of (in)equality constraints
    type :: CGStructureUDT

        ! Description
        !============
        ! UDT that contains information (for one constraint) of
        ! which dof groups (in 'dofgroups') can be used to which the
        ! (in)equality constraint can be attributed to.

        integer(I8), allocatable    :: dofgroups(:)

    end type

    ! Abstract types
    !===============
    ! Generic constraint type
    type, abstract, extends(ConstraintsUDT) :: GenericConstraintsGDUDT

        ! Description
        !============
        ! Generic type for grid deformation constraints. Inherits from
        ! the mother constraint type defined in gdmod_constraints.

        ! Fields:
        integer(I8)                 :: ncon = 0 ! number of constraints

    contains

        ! Initialization
        procedure(InitializeConstraintsINT), deferred :: Initialize

        ! Evaluation
        procedure(EvaluateConstraintsINT), deferred :: Evaluate

    end type

    ! Specific constraint types
    !==========================
    ! Flux function constraints
    type, extends(GenericConstraintsGDUDT) :: FluxfunctionConstraintsUDT

        ! Description
        !============
        ! Flux function constraints. Fixes the flux function values
        ! for a set of desired nodes. The following fields are added:
        ! - fixedpoints:    specifies which non-flux surface points
        !                   should have fixed psi value
        ! - specialpoints:  specifies which points should have a fixed
        !                   flux value that is equal to another point.
        !                   the first point in the ID field of this
        !                   structure is the flux value that is used
        !                   to fix the other points to. Useful for
        !                   separatrices etc of which the flux value
        !                   is a priori unknown (e.g. due to 'grad'
        !                   type of X-point constraints)
        ! - fluxsurfaces    specifies which vertices belong to a flux
        !                   surface. This is the classical alignment
        !                   constraint.
        ! - tangencypoints  specifies which points should be treated as
        !                   tangency points. These get an additional
        !                   constraint, namely that the normal on the
        !                   vessel boundary should be perpendicular to
        !                   the magnetic field vector.
        ! - ncon:           (inherited) total number of constraints

        ! No other routines than the standard initialization, evaluation
        ! and destruction routines are implemented nor needed.

        ! Fields:
        type(FFCStructureUDT), allocatable  :: fixedpoints(:)
        type(FFCStructureUDT), allocatable  :: specialpoints(:)
        type(FFCStructureUDT), allocatable  :: fluxsurfaces(:)
        type(FFCStructureUDT), allocatable  :: tangencypoints(:)
        integer(I8)                         :: nfixedpoints, &
            nspecialpoints, nfluxsurfaces, ntangencypoints, nspcon, &
            nfscon, ntpcon, nfpcon

    contains

        ! Initialization
        procedure :: Initialize     => InitializeFluxfunctionConstraints

        ! Evaluation
        procedure :: Evaluate       => EvaluateFluxfunctionConstraints
        procedure :: EvaluateDerivativesCoordinates     => &
            EvaluateCoordinatesDerivativesFluxFunctionConstraints
        procedure :: EvaluateDerivativesFlux            => &
            EvaluateFluxDerivativesFluxFunctionConstraints

        ! Data
        procedure :: WriteData      => WriteDataFluxfunctionConstraints

        ! Destructor
        final :: DestroyFluxfunctionConstraints

    end type


    ! Boundary function constraints
    type, extends(GenericConstraintsGDUDT) :: BoundaryFunctionConstraintsUDT

        ! Description
        !============
        ! Boundary function constraints. The boundary function is given
        ! by the 'psf' polygon (here, only one polygon can be specified
        ! though this may be extended in the future). The vertices on
        ! which these conditions have to be imposed are given in the
        ! 'vertID' array (number of vertices given in nvertID)

        ! Fields:
        class(PolygonLevelsetFunction2DUDT), allocatable    :: plf
        integer(I8), allocatable            :: vert(:)
        integer(I8)                         :: nvert

    contains

        ! Initialization
        procedure :: Initialize => InitializeBoundaryFunctionConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateBoundaryFunctionConstraints
        procedure :: EvaluateDerivativesCoordinates     => &
            EvaluateCoordinatesDerivativesBoundaryFunctionConstraints
        procedure :: EvaluateDerivativesFlux            => &
            EvaluateFluxDerivativesBoundaryFunctionConstraints
        procedure :: EvaluateVesselcoordinatesDerivativesBoundaryFunctionConstraints

        ! Update
        procedure :: Update     => UpdateBoundaryFunctionConstraints

    end type

    ! X-point constraints
    type, extends(GenericConstraintsGDUDT) :: XPointConstraintsUDT

        ! Description
        !============
        ! X-point constraints. Constrains the location of the x-point by
        ! fixing the initial x-point coordinates or look at where the
        ! gradient of the flux function vanishes.
        ! The following fields are added:
        !
        ! - xpind:      the x-point vertex IDs
        ! - nxpind:     the total number of x-points
        ! - locx/y:     x and y coordinate locations of the x-points
        ! - meth:       method to contrain x-point. Can be 'loc' or
        !               'grad'. In case of 'loc', the location is fixed
        !               based on locx/locy. Otherwise, gradient = 0
        !               conditions are imposed.

        ! Note 1: the initial x-point location and x-point vertices
        ! are determined by the initial grid. The x-point indices are
        ! determined using the DetermineXPoints routine in
        ! gdmod_utility_optimization.

        ! Fields:
        integer(I8), allocatable            :: xpind(:)
        integer(I8)                         :: nxpind
        real(R8), allocatable               :: locx(:)
        real(R8), allocatable               :: locy(:)
        character(:), allocatable           :: meth

    contains

        ! Initialization
        procedure :: Initialize => InitializeXPointConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateXPointConstraints
        procedure :: EvaluateCoordinatesDerivative &
            => EvaluateCoordinatesDerivativeXPointConstraints

    end type

    ! Edge lengths constraints
    type, extends(GenericConstraintsGDUDT) :: EdgeLengthsConstraintsUDT

        ! Description
        !============
        ! Edge length constraints. Constrain an arbitrary set of edges
        ! to have a certain edge length (in meters) - see also
        ! InitializeEdgeLengthsConstraints for additional predefined
        ! edge sets. The following fields are added:
        !
        ! - nedges:     (scalar) total number of edges to constrain
        ! - edgevert:   (nedges-by-2) edge vertex IDs
        ! - d           (nedges-by-1) desired length for each edge

        ! Fields:
        integer(I8), allocatable            :: edgevert(:, :)
        integer(I8)                         :: nedges
        real(R8), allocatable               :: d(:)

    contains

        ! Initialization
        procedure :: Initialize => InitializeEdgeLengthsConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateEdgeLengthsConstraints

    end type

    ! Orthogonality constraints
    type, extends(GenericConstraintsGDUDT) :: OrthogonalityConstraintsUDT

        ! Description
        !============
        ! Orthogonality constraints. These constraints fix certain edges
        ! to be orthogonal to the magnetic field. See the initialization
        ! routine (InitializeOrthogonalityConstraints) to see the
        ! different options that are available to determine which
        ! edges to constrain. Typically, edges near the core and SOL
        ! are constrained.
        ! The following fields are added:
        !
        ! - nedges:     (scalar) total number of edges to constrain
        ! - edgevert:   (nedges-by-2) edge vertex IDs
        ! - radiallines: polygon set with radial lines built from the
        !               edges

        ! Fields:
        integer(I8), allocatable            :: edgevert(:, :)
        integer(I8)                         :: nedges
        type(PolygonSetUDT)                 :: radiallines

    contains

        ! Initialization
        procedure :: Initialize => InitializeOrthogonalityConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateOrthogonalityConstraints

    end type

    ! Fixed flux value constraints
    type, extends(GenericConstraintsGDUDT)  :: FixedFluxvaluesConstraintsUDT

        ! Description
        !============
        ! Constraints to fix flux values explicitly. Only applicable
        ! if the design variables include the flux values. The following
        ! fields are specified:
        !
        ! - fsind       : indices of which flux surfaces are constrained
        ! - psiind      : indices of corresponding psi value design
        !               variable for each constraint
        ! - psid        : desired flux value of these flux surfaces
        !
        ! The options (see gdmod_userinput) allow to specify either
        ! automatically or manually the flux values at which core and
        ! outer flux surface(s) should be fixed.

        integer(I8), allocatable    :: fsind(:), psiind(:)
        real(R8), allocatable       :: psid(:)

    contains

        ! Initialization
        procedure :: Initialize => InitializeFixedFluxvaluesConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateFixedFluxvaluesConstraints


    end type

    ! Linefolding constraints
    type, extends(GenericConstraintsGDUDT)  :: LinefoldingConstraintsUDT

        ! Description
        !============
        ! Constraints to prevent overlapping coordinate lines and cells.
        ! There are three 'types' that can be used: 'poloidal', 'radial'
        ! and 'vessel'. Each has a different coordinate direction along
        ! which folding may not occur. This coordinate direction is
        ! given by the interpolant representation of the magnetic field
        ! and vessel levelset. Since line folding is prevented by
        ! imposing conditions to the inner product of the face tangent
        ! and coordinate vector, a small number has to be added to
        ! prevent coinciding vertices etc. This number should be
        ! sufficiently small (e.g. three orders of magnitude smaller
        ! than the smallest feature in the grid or something).
        ! Additionally, the constraints are ill-defined if both
        ! derivatives of the polygon levelset function are zero (e.g. at
        ! magnetic field extrema). To hedge for this, there is a field
        ! 'fieldtol' that specifies the absolute tolerance on the field
        ! vector under which the constraint is deactivated.

        ! Options basically
        real(R8)                    :: smallnumber, fieldtol

        ! Fields
        integer(I8)                 :: nvpairspol, nvpairsrad, &
            nvpairsves
        integer(I8), allocatable    :: vpairspol(:, :), vpairsrad(:, :), &
            vpairsves(:, :)
        real(R8), allocatable       :: signvecpol(:), signvecrad(:), &
            signvecves(:)

    contains

        ! Initialization
        procedure :: Initialize => InitializeLinefoldingConstraints

        ! Evaluation
        procedure :: Evaluate   => EvaluateLinefoldingConstraints

        ! Update
        procedure :: Update     => UpdateLinefoldingConstraints


    end type

    ! In-vessel constraints
    type, extends(BoundaryFunctionConstraintsUDT) :: InVesselConstraintsUDT

        ! Description
        !============
        ! (inequality) constraints to keep vertices inside the vessel
        ! structure. May help for sudden geometry changes etc. Should
        ! eventually boil down in combination with flux function
        ! constraints to moving the vertices (approximately) along the
        ! magnetic field lines. It is important to note that we assume
        ! that the vessel is closed and that the interior is defined
        ! by negative values of the vessel boundary representation.
        ! Therefore, we can inherit from the boundary constraints and
        ! we only need to overwrite the initialization of the
        ! constraints.

        ! No additional fields to define

    contains

        procedure :: Initialize         => InitializeInVesselConstraints

    end type

    ! Overarching types
    !==================
    ! Equality constraints
    type, extends(ConstraintsUDT) :: EqConGDUDT

        ! Description
        !============
        ! This type contains all the different constraints as different
        ! derived types. For each type of constraint, a different
        ! type is defined.

        ! Total number of constraints
        integer(I8)                         :: neqcon = 0

        ! Constraint switches
        logical                             :: dofluxfunction = .false.
        logical                             :: doboundaryfunction = .false.
        logical                             :: doxpoints = .false.
        logical                             :: doedgelengths = .false.
        logical                             :: doorthogonality = .false.
        logical                             :: dofixedfluxvalues = .false.

        type(FluxfunctionConstraintsUDT)    :: fluxfunction
        type(BoundaryFunctionConstraintsUDT):: boundaryfunction
        type(XPointConstraintsUDT)          :: xpoints
        type(EdgeLengthsConstraintsUDT)     :: edgelengths
        type(OrthogonalityConstraintsUDT)   :: orthogonality
        type(FixedFluxvaluesConstraintsUDT) :: fixedfluxvalues

    contains

        ! Procedure to initialize constraints
        procedure :: Initialize         => InitializeEqCon

        ! Procedure to evaluate constraints
        procedure :: Evaluate           => EvaluateEqCon

    end type

    ! Inequality constraints
    type, extends(ConstraintsUDT) :: IneqConGDUDT

        ! Total number of inequality constraints
        integer(I8)             :: nineqcon = 0

        ! Constraint switches
        logical                             :: dolinefolding = .false.
        logical                             :: doinvessel = .false.

        type(LinefoldingConstraintsUDT)     :: linefolding
        type(InVesselConstraintsUDT)        :: invessel

    contains

        ! Constraints initialization
        procedure :: Initialize         => InitializeIneqCon

        ! Constraints evaluation
        procedure :: Evaluate           => EvaluateIneqCon

    end type

    ! All constraints for the grid deformation
    type, extends(ConstraintsUDT) :: ConstraintsGDUDT

        ! Description
        !============
        ! Defines the basic optimization problem: it has a set of
        ! design variables, constraints, and a cost function.

        ! Fields:

        ! Equality constraints
        type(EqConGDUDT)        :: eqcon

        ! Inequality constraints
        type(IneqConGDUDT)      :: ineqcon

    contains

        ! Initialization
        procedure :: Initialize         => InitializeConstraints

        ! Number of constraints getter
        procedure :: GetConstraintsDimensions

        ! Evaluation
        ! procedure :: Evaluate           => EvaluateEqualityConstraintsGD

        ! Housekeeping

    end type

    ! Monitor
    !========
    type ConstraintsMonitorUDT

        ! Description
        !============
        ! This object can be used to keep track of the amount of
        ! constraints that are imposed per vertex. As such, it can be
        ! checked whether the problem will be overly constrained and to
        ! prevent constraints to be imposed on certain vertices if
        ! a maximum number is reached. For now, the equality constraints
        ! per vertex are counted in eqvcc, the inequality in ineqvcc.
        ! The maximum number of equality and inequality constraints per
        ! vertex are set in maxeqvcc and maxineqvcc. The methods are
        ! the following:

        ! - Initialize: initialization routine that takes the grid,
        !   magnetic field, and environment as input (in case of
        !   extension of this object in the future)
        ! - CheckConstraints: will check if the problem is already
        !   overly constrained, or may become overly constrained if
        !   inequality constraints become active.

        ! Fields:
        integer(I8), allocatable           :: eqvcc(:), ineqvcc(:)
        integer(I8), allocatable           :: maxeqvcc(:), maxineqvcc(:)

    contains

        ! Routines
        procedure :: Initialize         => InitializeMonitor
        ! procedure :: CheckConstraints   => CheckConstraintsMonitor

        ! Finalization
        final :: DestroyMonitor

    end type

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    abstract interface

        ! Constraint initialization
        subroutine InitializeConstraintsINT(constraints, grid, &
            magneticField, environment, monitor, designvariables, &
            options)

            ! Description
            !============
            ! This routine serves as general initialization routine for
            ! a generic grid deformation constraint (that inherits from
            ! the GenericConstraintUDT type)

            ! Import
            import :: GenericConstraintsGDUDT, GridUDT, &
                MagneticFieldUDT, EnvironmentUDT, &
                ConstraintsMonitorUDT, ConstraintOptionsUDT, &
                DesignVariablesGDUDT

            ! Declare
            class(GenericConstraintsGDUDT)      :: constraints
            type(GridUDT)                       :: grid
            type(MagneticFieldUDT)              :: magneticField
            type(EnvironmentUDT)                :: environment
            type(ConstraintsMonitorUDT)         :: monitor
            type(ConstraintOptionsUDT)          :: options
            class(DesignVariablesGDUDT)         :: designvariables

        end subroutine

        ! Constraint evaluation
        subroutine EvaluateConstraintsINT(constraints, G, gradG, &
            hessG, grid, magneticField, environment, &
            dogradient, dohessian, designvariables, lambda, varin, &
            valuesin, dGdvarin, dgradGdvarin)

            ! Description
            !============
            ! This reoutine serves as a general evaluation routine for
            ! a generic grid deformation constraint.

            ! Import
            import :: GenericConstraintsGDUDT, MySparseUDT, GridUDT, &
                R8, MagneticFieldUDT, EnvironmentUDT, &
                DesignVariablesGDUDT

            ! Declare
            class(GenericConstraintsGDUDT)  :: constraints
            real(R8), allocatable           :: G(:), lambda(:)
            type(MySparseUDT)               :: hessG, gradG
            type(GridUDT)                   :: grid
            type(MagneticFieldUDT)          :: magneticField
            type(EnvironmentUDT)            :: environment
            logical                         :: dogradient, dohessian
            class(DesignVariablesGDUDT)     :: designvariables

            character(*), intent(in), optional  :: varin
            real(R8), intent(in), optional      :: valuesin(:)
            type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                           CONSTRAINTS MONITOR                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeMonitor(monitor, grid, magneticField, &
        environment, designvariables)

        ! Description
        !============
        ! Initializes the constraints monitor structure. It is assumed
        ! that the grid, magnetic field and environment are properly
        ! allocated and initialized.

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsMonitorUDT)        :: monitor
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        class(DesignVariablesGDUDT)         :: designvariables

        ! Loop variables

        ! Auxiliary

        ! Initialize
        !===========
        ! Allocate
        allocate(monitor%eqvcc(grid%vert%ntot), monitor%maxeqvcc(grid%vert%ntot))
        allocate(monitor%ineqvcc(grid%vert%ntot), monitor%maxineqvcc(grid%vert%ntot))

        ! Initialize
        monitor%eqvcc(:)        = 0
        monitor%ineqvcc(:)      = 0
        monitor%maxeqvcc(:)     = 2 ! for most cases this is fine
        monitor%maxineqvcc(:)   = 1000 ! a stupid large number - can impose any number

        ! Check
        select type (designvariables)

        type is (DesignVariablesCoordinatesFluxUDT)

            monitor%maxeqvcc(:) = 3

        class default

        end select

    end subroutine

    ! Destructor
    subroutine DestroyMonitor(monitor)

        ! Description
        !============
        ! Finalization for the constraints monitor object

        ! Declare variables
        !==================
        ! Arguments
        type(ConstraintsMonitorUDT)         :: monitor

        ! Destroy
        !========
        deallocate(monitor%eqvcc, monitor%ineqvcc)

    end subroutine

    !------------------------------------------------------------------!
    !                           GENERAL CONSTRAINTS                    !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeConstraints(constraints, grid, magneticField, &
        environment, designvariables, options)

        ! Description
        !============
        ! Routine that initializes the equality and inequality
        ! constraints, using the initialization routines of those
        ! objects.

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)     :: constraints
        type(GridUDT)               :: grid
        type(MagneticFieldUDT)      :: magneticField
        type(EnvironmentUDT)        :: environment
        type(ConstraintOptionsUDT)  :: options
        type(ConstraintsMonitorUDT) :: monitor
        class(DesignVariablesGDUDT) :: designvariables

        ! Loop variables

        ! Auxiliary variables

        ! Initialize monitor
        !===================
        call monitor%Initialize(grid, magneticField, environment, designvariables)

        ! Initialize constraints
        !=======================
        ! Equality constraints
        call constraints%eqcon%Initialize(grid, magneticField, &
            environment, options, designvariables, monitor)

        ! Inequality constraints
        call constraints%ineqcon%Initialize(grid, magneticField, &
            environment, options, designvariables, monitor)

    end subroutine

    ! Dimension getter
    subroutine GetConstraintsDimensions(constraints, neqcon, nineqcon)

        ! Description
        !============
        ! Return the current dimensions of the equality and inequality
        ! constraints. Can be used for initialization of other
        ! quantities at higher levels.

        ! Declare variables
        !==================
        ! Arguments
        class(ConstraintsGDUDT)         :: constraints
        integer(I8), intent(out)        :: neqcon, nineqcon

        ! Extract dimensions
        !===================
        ! Stored in eqcon, ineqcon
        neqcon = constraints%eqcon%neqcon
        nineqcon = constraints%ineqcon%nineqcon

    end subroutine



    !------------------------------------------------------------------!
    !                           EQUALITY CONSTRAINTS                   !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeEqCon(constraints, grid, magneticField, &
        environment, constraintoptions, designvariables, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)           :: constraints
        type(GridUDT)               :: grid
        type(MagneticFieldUDT)      :: magneticField
        type(EnvironmentUDT)        :: environment
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor
        class(DesignVariablesGDUDT) :: designvariables

        ! Loop variables

        ! Auxiliary variables

        ! Initialize constraints
        !=======================
        constraints%neqcon = 0

        ! X-points
        if (constraintoptions%xpoints == 1) then

            ! Set the logical
            constraints%doxpoints = .true.

            ! Initialize
            call constraints%xpoints%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%xpoints%ncon

            ! Print
            print *, 'number of x-point function constraints: ', &
                constraints%xpoints%ncon


        else
            ! Set to false, don't initialize
            constraints%doxpoints = .false.

        end if

        ! Boundary function
        if (constraintoptions%boundaryfunctions == 1) then

            ! Set the logical
            constraints%doboundaryfunction = .true.

            ! Initialize
            call constraints%boundaryfunction%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%boundaryfunction%ncon

            ! Print
            print *, 'number of boundary function constraints: ', &
                constraints%boundaryfunction%ncon


        else
            ! Set to false, don't initialize
            constraints%doboundaryfunction = .false.

        end if

        ! Flux function
        if (constraintoptions%fluxfunction == 1) then
            ! Set the logical
            constraints%dofluxfunction = .true.

            ! Initialize
            call constraints%fluxfunction%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fluxfunction%ncon

            ! Print
            print *, 'number of flux function constraints: ', &
                constraints%fluxfunction%ncon

        else
            ! Set to false, don't initialize
            constraints%dofluxfunction = .false.

        end if

        ! Edge lengths
        if (constraintoptions%edgelengths == 1) then
            ! Set the logical
            constraints%doedgelengths = .true.

            ! Initialize
            call constraints%edgelengths%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%edgelengths%ncon

            ! Print
            print *, 'number of edge lengths constraints: ', &
                constraints%edgelengths%ncon

        else
            ! Set to false, don't initialize
            constraints%doedgelengths = .false.

        end if

        ! Orthogonality
        if (constraintoptions%orthogonality == 1) then
            ! Set the logical
            constraints%doorthogonality = .true.

            ! Initialize
            call constraints%orthogonality%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%orthogonality%ncon

            ! Print
            print *, 'number of orthogonality constraints: ', &
                constraints%orthogonality%ncon

        else
            ! Set to false, don't initialize
            constraints%doorthogonality = .false.

        end if

        ! Fixed flux values
        if (constraintoptions%fixedfluxvalues == 1) then
            ! Set the logical
            constraints%dofixedfluxvalues = .true.

            ! Initialize
            call constraints%fixedfluxvalues%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%neqcon = constraints%neqcon + &
                constraints%fixedfluxvalues%ncon

            ! Print
            print *, 'number of fixed flux value constraints: ', &
                constraints%fixedfluxvalues%ncon

        else
            ! Set to false, don't initialize
            constraints%dofixedfluxvalues = .false.

        end if



    end subroutine

    ! Constraint evaluation
    subroutine EvaluateEqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables, lambda, varin, valuesin, dGdvarin, dgradGdvarin)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is
        ! checked whether it is imposed, and the contributions are
        ! added by calling the evaluation routine of each constraint.

        ! Declare variables
        !==================
        ! Arguments
        class(EqConGDUDT)               :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian
        class(DesignVariablesGDUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop
        integer(I8)                     :: ic, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary
        real(R8), allocatable           :: G_flux(:), lambda_flux(:)
        type(MySparseUDT)               :: gradG_flux, hessG_flux, &
            dG_fluxdvar, dgradG_fluxdvar

        real(R8), allocatable           :: G_bnd(:), lambda_bnd(:)
        type(MySparseUDT)               :: gradG_bnd, hessG_bnd, &
            dG_bnddvar, dgradG_bnddvar

        real(R8), allocatable           :: G_xp(:), lambda_xp(:)
        type(MySparseUDT)               :: gradG_xp, hessG_xp, &
            dG_xpdvar, dgradG_xpdvar

        real(R8), allocatable           :: G_el(:), lambda_el(:)
        type(MySparseUDT)               :: gradG_el, hessG_el, &
            dG_eldvar, dgradG_eldvar

        real(R8), allocatable           :: G_orth(:), lambda_orth(:)
        type(MySparseUDT)               :: gradG_orth, hessG_orth, &
            dG_orthdvar, dgradG_orthdvar

        real(R8), allocatable           :: G_ffv(:), lambda_ffv(:)
        type(MySparseUDT)               :: gradG_ffv, hessG_ffv, &
            dG_ffvdvar, dgradG_ffvdvar

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if

        ! Set the constraint counter
        ic = 0

        ! Constraint values
        !==================
        gradG = SpZeros(designvariables%nphi, 0) ! will concatenate
        hessG = SpZeros(designvariables%nphi, designvariables%nphi) ! will add
        dGdvar = SpZeros(0, size(values, 1))
        dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        ! X-point constraints
        !--------------------
        if (constraints%doxpoints) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%xpoints%ncon)]

            ! Allocate & initialize
            lambda_xp = lambda(conindex)

            ! Call the evaluation routine
            call constraints%xpoints%Evaluate(G_xp, &
                gradG_xp, hessG_xp, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_xp, var, values, dG_xpdvar, dgradG_xpdvar)

            ! Assign
            G(conindex) = G_xp
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_xp, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_xp
            end if
            dGdvar = dGdvar%Concatenate(dG_xpdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_xpdvar

            ! Update the constraint counter
            ic = ic + constraints%xpoints%ncon

        end if

        ! Boundary function constraints
        !------------------------------
        if (constraints%doboundaryfunction) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%boundaryfunction%ncon)]

            ! Allocate & initialize
            lambda_bnd = lambda(conindex)

            ! Call the evaluation routine
            call constraints%boundaryfunction%Evaluate(G_bnd, &
                gradG_bnd, hessG_bnd, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_bnd, var, values, dG_bnddvar, dgradG_bnddvar)

            ! Assign
            G(conindex) = G_bnd
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_bnd, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_bnd
            end if
            dGdvar = dGdvar%Concatenate(dG_bnddvar, 1)
            dgradGdvar = dgradGdvar + dgradG_bnddvar

            ! Update the constraint counter
            ic = ic + constraints%boundaryfunction%ncon

        end if


        ! Flux function constraints
        !--------------------------
        if (constraints%dofluxfunction) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%fluxfunction%ncon)]

            ! Allocate & initialize
            lambda_flux = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fluxfunction%Evaluate(G_flux, &
                gradG_flux, hessG_flux, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_flux, var, values, dG_fluxdvar, dgradG_fluxdvar)

            ! Assign
            G(conindex) = G_flux
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_flux, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_flux
            end if
            dGdvar = dGdvar%Concatenate(dG_fluxdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_fluxdvar

            ! Update the constraint counter
            ic = ic + constraints%fluxfunction%ncon

        end if

        ! Edge lengths constraints
        !-------------------------
        if (constraints%doedgelengths) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%edgelengths%ncon)]

            ! Allocate & initialize
            lambda_el = lambda(conindex)

            ! Call the evaluation routine
            call constraints%edgelengths%Evaluate(G_el, &
                gradG_el, hessG_el, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_el, var, values, dG_eldvar, dgradG_eldvar)

            ! Assign
            G(conindex) = G_el

            ! Update the gradient column indices
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_el, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_el
            end if
            dGdvar = dGdvar%Concatenate(dG_eldvar, 1)
            dgradGdvar = dgradGdvar + dgradG_eldvar

            ! Update the constraint counter
            ic = ic + constraints%edgelengths%ncon

        end if

        ! Orthogonality constraints
        !--------------------------
        if (constraints%doorthogonality) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%orthogonality%ncon)]

            ! Allocate & initialize
            lambda_orth = lambda(conindex)

            ! Call the evaluation routine
            call constraints%orthogonality%Evaluate(G_orth, &
                gradG_orth, hessG_orth, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_orth, var, values, dG_orthdvar, dgradG_orthdvar)

            ! Assign
            G(conindex) = G_orth
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_orth, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_orth
            end if
            dGdvar = dGdvar%Concatenate(dG_orthdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_orthdvar

            ! Update the constraint counter
            ic = ic + constraints%orthogonality%ncon

        end if

        ! Fixed flux values constraints
        !------------------------------
        if (constraints%dofixedfluxvalues) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%fixedfluxvalues%ncon)]

            ! Allocate & initialize
            lambda_ffv = lambda(conindex)

            ! Call the evaluation routine
            call constraints%fixedfluxvalues%Evaluate(G_ffv, &
                gradG_ffv, hessG_ffv, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_ffv, var, values, dG_ffvdvar, dgradG_ffvdvar)

            ! Assign
            G(conindex) = G_ffv
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_ffv, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_ffv
            end if
            dGdvar = dGdvar%Concatenate(dG_ffvdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_ffvdvar

            ! Update the constraint counter
            ic = ic + constraints%fixedfluxvalues%ncon

        end if

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                          INEQUALITY CONSTRAINTS                  !
    !------------------------------------------------------------------!
    ! Initialization
    subroutine InitializeIneqCon(constraints, grid, magneticField, &
        environment, constraintoptions, designvariables, monitor)

        ! Description
        !============
        ! Routine that initializes the desired constraints, based on the
        ! user options d

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)         :: constraints
        type(GridUDT)               :: grid
        type(MagneticFieldUDT)      :: magneticField
        type(EnvironmentUDT)        :: environment
        type(ConstraintOptionsUDT)  :: constraintoptions
        type(ConstraintsMonitorUDT) :: monitor
        class(DesignVariablesGDUDT)  :: designvariables

        ! Loop variables

        ! Auxiliary variables

        ! Initialize
        !===========

        ! Initialize constraints
        !=======================
        constraints%nineqcon = 0

        ! Linefolding
        if (constraintoptions%linefolding == 1) then

            ! Set the logical
            constraints%dolinefolding = .true.

            ! Initialize
            call constraints%linefolding%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%nineqcon = constraints%nineqcon + &
                constraints%linefolding%ncon

            ! Print
            print *, 'number of linefolding constraints: ', &
                constraints%linefolding%ncon


        else
            ! Set to false, don't initialize
            constraints%dolinefolding = .false.

        end if

        ! Invessel
        if (constraintoptions%invessel == 1) then

            ! Set the logical
            constraints%doinvessel = .true.

            ! Initialize
            call constraints%invessel%Initialize(grid, &
                magneticField, environment, monitor, designvariables, &
                constraintoptions)

            ! Add constraints number
            constraints%nineqcon = constraints%nineqcon + &
                constraints%invessel%ncon

            ! Print
            print *, 'number of invessel constraints: ', &
                constraints%invessel%ncon


        else
            ! Set to false, don't initialize
            constraints%doinvessel = .false.

        end if

    end subroutine

    ! Constraint evaluation
    subroutine EvaluateIneqCon(constraints, G, gradG, hessG, &
        grid, magneticField, environment, dogradient, dohessian, &
        designvariables, lambda, varin, valuesin, dGdvarin, dgradGdvarin)

        ! Description
        !============
        ! This routine evaluates the constraints G and the corresponding
        ! gradient and hessian. To do so, every type of constraint is
        ! checked whether it is imposed, and the contributions are
        ! added by calling the evaluation routine of each constraint.

        ! Declare variables
        !==================
        ! Arguments
        class(IneqConGDUDT)             :: constraints
        real(R8), intent(inout)         :: G(:)
        real(R8), intent(in)            :: lambda(:)
        type(MySparseUDT)               :: gradG, hessG
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField
        type(EnvironmentUDT)            :: environment
        logical                         :: dogradient, dohessian
        class(DesignVariablesGDUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop
        integer(I8)                     :: ic, k
        integer(I8), allocatable        :: conindex(:)

        ! Auxiliary
        real(R8), allocatable           :: G_lf(:), lambda_lf(:)
        type(MySparseUDT)               :: gradG_lf, hessG_lf, &
            dG_lfdvar, dgradG_lfdvar

        real(R8), allocatable           :: G_iv(:), lambda_iv(:)
        type(MySparseUDT)               :: gradG_iv, hessG_iv, &
            dG_ivdvar, dgradG_ivdvar

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if

        ! Set the constraint counter
        ic = 0

        ! Set derivatives

        ! Constraint values
        !==================
        gradG = SpZeros(designvariables%nphi, 0) ! will concatenate
        hessG = SpZeros(designvariables%nphi, designvariables%nphi) ! will add
        dGdvar = SpZeros(0, size(values, 1))
        dgradGdvar = SpZeros(designvariables%nphi, size(values))


        ! Invessel constraints
        !---------------------
        if (constraints%doinvessel) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%invessel%ncon)]

            ! Allocate & initialize
            lambda_iv = lambda(conindex)

            ! Call the evaluation routine
            call constraints%invessel%Evaluate(G_iv, &
                gradG_iv, hessG_iv, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_iv, var, values, dG_ivdvar, dgradG_ivdvar)

            ! Assign
            G(conindex) = G_iv
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_iv, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_iv
            end if
            dGdvar = dGdvar%Concatenate(dG_ivdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_ivdvar

            ! Update the constraint counter
            ic = ic + constraints%invessel%ncon

        end if

        ! Linefolding constraints
        !------------------------
        if (constraints%dolinefolding) then
            ! Construct the constraint index
            conindex = [(k, k = ic+1, ic+constraints%linefolding%ncon)]

            ! Allocate & initialize
            lambda_lf = lambda(conindex)

            ! Call the evaluation routine
            call constraints%linefolding%Evaluate(G_lf, &
                gradG_lf, hessG_lf, &
                grid, magneticField, environment, dogradient, &
                dohessian, designvariables, &
                lambda_lf, var, values, dG_lfdvar, dgradG_lfdvar)

            ! Assign
            G(conindex) = G_lf
            if (dogradient) then
                gradG = gradG%Concatenate(gradG_lf, 2)
            end if
            if (dohessian) then
                hessG = hessG + hessG_lf
            end if
            dGdvar = dGdvar%Concatenate(dG_lfdvar, 1)
            dgradGdvar = dgradGdvar + dgradG_lfdvar

            ! Update the constraint counter
            ic = ic + constraints%linefolding%ncon

        end if

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if


    end subroutine


    !------------------------------------------------------------------!
    !                           FLUX FUNCTION                          !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeFluxfunctionConstraints(constraints, grid, &
        magneticField, environment, monitor, designvariables, options)

        ! Description
        !============
        ! Initialize the required fields of the flux function
        ! constraints. These constraints impose for each vertex that
        ! lies on a flux surface (i.e. has a flux surface ID which is
        ! non-zero) by fixing its flux value. The flux values can be
        ! determined in different ways. Here, we compute the initial
        ! psi values by averaging the current psi values of the vertices
        ! that lie on a flux surface. This hedges a bit for
        ! discretization errors originating from the re-interpretation
        ! of the magnetic field data as an interpolant here
        ! instead of a linear interpolant in the grid generator.

        ! There are several different 'flavours' on how these
        ! constraints can be imposed. For classical flux surfaces, all
        ! data is stored in the structure 'fluxsurfaces'. Separatrices
        ! are dealt with using the 'special points', where the first
        ! point is the X-point, and the flux values of the other points
        ! are equated to this one. 'tangencypoints' specifies which
        ! points have a magnetic field vector that should be (at the
        ! solution) tangential to the vessel boundary vector. One should
        ! be careful with this constraint - a solution may not exist!
        ! Different treatment options are therefore available through
        ! the option 'tangencypointtreatment':
        !   - 'tangencypoint': classic tangency point constraint.
        !   If the design variables are only coordinates, then the
        !   tangency points are considered as special points, since
        !   the additional constraint fixes its position, similar to
        !   X-points. If not, we can additionally impose aligment.
        !   Combined with the boundary constraints, this triplet fixes
        !   the coordinates and flux values
        !   - 'align': treat as a vertex of a classical flux surface.
        !   One should be careful as this may lead to constraint
        !   singularity (at an intermediate point or solution, the
        !   vessel boundary constraints and flux value constraint may be
        !   collinear)
        !   - 'noconstraint': no additional flux surface constraint is
        !   imposed. This may lead to local non-alignment of the
        !   adjacent faces, but may stabilize the problem.

        ! Notes
        !======
        ! Note 1: the last argument of this function is a derived type
        ! used to monitor the constraints. It can be used to make sure
        ! that the problem is - at least not by the equality constraints
        ! - is overly constrained. It is UP TO THE DEVELOPER to use this
        ! monitor properly!

        ! Note 2: the constraint options are passed to this function in
        ! order to determine whether e.g. boundary nodes should be
        ! considered for the constraints.

        ! Note 3: for vessel mode grids (grids reaching up to the
        ! vessel wall), some vertices may not belong to a flux surface.
        ! To fix them anyway on the vessel wall, the logical
        ! 'fixfarvesselflux' can be set to 'true'. This is the
        ! recommended default value.

        ! Note 4: for target mode grids, the flux value at the 'corners'
        ! is typically also determined by averaging over the field line.
        ! However, this gives issues in some cases where the targets are
        ! nearly flux-aligned. To avoid this, set the logical
        ! 'fixfluxalignedtargets' to true. This is the recommended
        ! default value for narrow grids. For wide grids, this should be
        ! false.

        ! Note 5: it is assumed that the grid contains the flux surface
        ! data and that the flux surfaces are numbered from 1 to nFs

        ! Initialize
        !===========
        ! Modules
        use gdmod_plots
        use mod_definitions, only : targetID, vesselID
        implicit none

        ! Declare variables
        !==================
        ! Arguments
        class(FluxfunctionConstraintsUDT)       :: constraints
        type(GridUDT)                           :: grid
        type(MagneticFieldUDT)                  :: magneticField
        type(EnvironmentUDT)                    :: environment
        type(ConstraintsMonitorUDT)             :: monitor
        type(ConstraintOptionsUDT)              :: options
        class(DesignVariablesGDUDT)              :: designvariables

        ! Loop variables
        integer(I8)                 :: i, k

        ! Auxiliary variables
        real(R8)                    :: tpsi
        real(R8), allocatable       :: PsiD_tmp(:), temppsi(:), &
            boxx(:, :), boxy(:, :)

        integer(I8)                 :: nxpind, ntpind, nfsIDs, nspc
        integer(I8), allocatable    :: vert_tmp(:), vertID(:), &
            order(:), xpind(:), tptype(:), fsxpind(:), fstpind(:), &
            allIDs(:), fsIDs(:), vID(:), tpind(:), tv(:)

        logical                     :: fixxp, fixtp
        logical, allocatable        :: delind(:), mask(:), isdouble(:), &
            isxp(:), isconstrainedv(:), isvesselvertex(:), &
            isvesselface(:), isfixedpoint(:), includev(:), istp(:), &
            isfsxp(:), isfstp(:), isnofieldline(:), istempfree(:)

        ! Initialize
        !===========
        ! Check allocation status
        if (allocated(constraints%tangencypoints)) then
            deallocate(constraints%tangencypoints)
        end if
        if (allocated(constraints%specialpoints)) then
            deallocate(constraints%specialpoints)
        end if
        if (allocated(constraints%fixedpoints)) then
            deallocate(constraints%fixedpoints)
        end if
        if (allocated(constraints%fluxsurfaces)) then
            deallocate(constraints%fluxsurfaces)
        end if

        ! Number of constraints
        constraints%ncon = 0

        ! Associate
        associate(&
            vert    => grid%vert,       x       => grid%vert%x,     &
            y       => grid%vert%y,     cc      => monitor%eqvcc,   &
            nv      => grid%vert%ntot,                              &
            maxcc   => monitor%maxeqvcc,                            &
            opt     => options%ffoptions,                           &
            fieldlineID             => grid%vert%fieldlineID,       &
            fixfluxalignedtargets   => options%ffoptions%fixfluxalignedtargets, &
            fixfarvesselflux        => options%ffoptions%fixfarvesselflux,      &
            fixallvesselvertices    => options%ffoptions%fixallvesselvertices,  &
            doboxoverride           => options%ffoptions%doboxoverride,         &
            nsp                     => constraints%nspecialpoints,  &
            nspcon                  => constraints%nspcon,          &
            nfs                     => constraints%nfluxsurfaces,   &
            nfscon                  => constraints%nfscon,          &
            nfp                     => constraints%nfixedpoints,    &
            nfpcon                  => constraints%nfpcon,          &
            ntp                     => constraints%ntangencypoints, &
            ntpcon                  => constraints%ntpcon           &
            )

        ! Allocate
        allocate(vert_tmp(nv),PsiD_tmp(nv), &
            temppsi(nv), mask(nv), delind(nv), &
            vertID(nv), istp(nv), &
            isconstrainedv(nv), isfixedpoint(nv), isnofieldline(nv), &
            istempfree(nv))

        ! Initialize
        PsiD_tmp(:) = 0
        temppsi(:)  = 0
        mask(:)     = .false.
        vID = [(i, i = 1, nv)]
        vert_tmp    = vID
        nfsIDs = maxval(vert%fieldlineID)
        allIDs = [(i, i = 1, nfsIDs)]
        isconstrainedv(:) = .false.
        isfixedpoint(:) = .false.

        nsp = 0
        ntp = 0
        nfp = 0
        nfs = 0

        ! Evaluate
        call magneticField%interp%Evaluate(x, y, 0, 0, PsiD_tmp)

        ! Check options
        select case (designvariables%type)

        case ('coordinates_desiredflux')

            ! Don't explicitly fix x-point
            fixxp = .false.

        case default

            ! Fix x-point
            fixxp = .true.

        end select

        ! Extract x-point indices
        call DetermineXPoints(xpind, nxpind, order, grid)
        allocate(isxp(nv))
        isxp(:)     = .false.
        isxp(xpind) = .true.

        ! Extract tangency point indices
        call DetermineTangencyPoints(tpind, ntpind, tptype, grid)
        istp(:)     = .false.
        istp(tpind) = .true.

        ! Check how to treat tangency points
        fixtp = .false.
        select case (opt%tangencypointtreatment)

        case ('tangencypoint')

            ! Check how many constraints there are
            do i = 1, ntpind

                if (cc(tpind(i)) < maxcc(tpind(i))) then
                    ! Update counter
                    ntp = ntp + 1
                end if
            end do

            ! Impose true tangency point constraints
            allocate(constraints%tangencypoints(ntp))
            do i = 1, ntpind
                if (cc(tpind(i)) < maxcc(tpind(i))) then
                    ! Impose
                    allocate(constraints%tangencypoints(i)%ID(1))
                    constraints%tangencypoints(i)%ID = tpind(i)
                    constraints%tangencypoints(i)%nID = 1
                    constraints%tangencypoints(i)%fsID = fieldlineID(tpind(i))

                    ! Update constraint counter
                    cc(tpind(i)) = cc(tpind(i)) + 1
                end if
            end do

            ! Check design variables
            select case (designvariables%type)

            case ('coordinates_desiredflux')

                ! Also impose aligment constraints
                istp(:) = .false.

            case default

                ! Treat as special point
                fixtp = .true.

            end select

        case ('align')

            ! classic alignment, don't do anything and set istp to false
            ntp = 0
            allocate(constraints%tangencypoints(ntp))
            istp(:) = .false.

        case ('noconstraint')

            ! Do not constrain, not even with aligment -> istp is still true
            ntp = 0
            allocate(constraints%tangencypoints(ntp))

        case default

            ! Throw error
            call gdErrorHandler('InitializeConstraintParameters: unknown option for tangency points')

        end select

        ! Special points
        !===============
        ! Get flux surfaces - assumed each X-point has its own surface,
        ! so connected double null not yet supported

        ! Get IDs of x-point and tangency point flux surfaces
        fsxpind = fieldlineID(xpind)
        fstpind = fieldlineID(tpind)
        allocate(isfsxp(nfsIDs), isfstp(nfsIDs))
        isfsxp(:) = .false.
        isfstp(:) = .false.
        isfsxp(fsxpind) = .true.
        do i = 1, size(fstpind)
            if (fstpind(i) /= 0) then
                isfstp(fstpind(i)) = .true.
            end if
        end do

        ! Check for doubles in x-point flux surfaces
        allocate(isdouble(nxpind))
        isdouble(:) = .false.
        do i = 1, nxpind
            if (isdouble(i)) then
                call gdErrorHandler('InitializeConstraintParametersFluxFunction: ' //&
                    ' x-points with same flux surface ID detected, not supported')
            end if
            isdouble(i) = .true.
        end do
        deallocate(isdouble)

        ! Check for doubles in tangency point flux surfaces
        allocate(isdouble(ntpind))
        isdouble(:) = .false.
        do i = 1, ntpind
            if (isdouble(i)) then
                call gdErrorHandler('InitializeConstraintParametersFluxFunction: ' //&
                    ' tangency points with same flux surface ID detected, not supported')
            end if
            isdouble(i) = .true.
        end do
        deallocate(isdouble)

        ! Exclude separatrices and/or tangency point flux surfaces?
        if (fixxp .and. fixtp) then
            ! Exclude
            allocate(fsIDs(count(.not. (isfsxp .or. isfstp))))
            fsIDs = pack(allIDs, (.not. (isfsxp .or. isfstp)))
        elseif (fixxp) then
            allocate(fsIDs(count(.not. isfsxp)))
            fsIDs = pack(allIDs, (.not. isfsxp))
        elseif (fixtp) then
            allocate(fsIDs(count(.not. isfstp)))
            fsIDs = pack(allIDs, (.not. isfstp))
        else
            fsIDs = allIDs
        end if

        ! Determine total number of special points
        if (fixxp) then
            nsp = nsp + nxpind
        end if
        if (fixtp) then
            ! Only type one tangency points are special points!
            nsp = nsp + count(tptype == 1)
        end if

        ! Allocate
        allocate(constraints%specialpoints(nsp))
        ! Loop over all X-points
        nspc = 0
        if (fixxp) then
            ! Update number of special points
            do i = 1, nxpind

                ! Update counter
                nspc = nspc + 1

                ! Add all vertices on this flux surface (xp first)
                allocate(tv(count(grid%vert%fieldlineID == fsxpind(i))))
                tv = pack(vID, grid%vert%fieldlineID == fsxpind(i))
                constraints%specialpoints(nspc)%ID = [xpind(i), pack(tv, tv .ne. xpind(i))]
                constraints%specialpoints(nspc)%nID = size(constraints%specialpoints(nspc)%ID, 1)
                constraints%specialpoints(nspc)%fsID = fsxpind(i)
                isconstrainedv(tv) = .true.
                cc(tv) = cc(tv) + 1
                deallocate(tv)

            end do
        end if

        ! Loop over all tangency points
        if (fixtp) then
            do i = 1, ntpind
                if (tptype(i) == 1) then
                    ! Update counter
                    nspc = nspc + 1

                    ! Add all vertices on this flux surface (tp first)
                    if (fstpind(i) /= 0) then
                        allocate(tv(count(grid%vert%fieldlineID == fstpind(i))))
                        tv = pack(vID, grid%vert%fieldlineID == fstpind(i))
                        constraints%specialpoints(nspc)%ID = [tpind(i), pack(tv, tv .ne. tpind(i))]
                        constraints%specialpoints(nspc)%nID = size(constraints%specialpoints(nspc)%ID, 1)
                        constraints%specialpoints(nspc)%fsID = fstpind(i)
                        isconstrainedv(tv) = .true.
                        cc(tv) = cc(tv) + 1
                        deallocate(tv)
                    else
                        tv = tpind(i:i)
                        constraints%specialpoints(nspc)%ID = [tv]
                        constraints%specialpoints(nspc)%nID = size(constraints%specialpoints(nspc)%ID, 1)
                        constraints%specialpoints(nspc)%fsID = fstpind(i)
                        isconstrainedv(tv) = .true.
                        cc(tv) = cc(tv) + 1
                        deallocate(tv)
                    end if
                end if
            end do
        end if

        ! Fixed points
        !=============
        ! Initialize
        mask(:) = .false.
        isfixedpoint(:) = .false.

        ! Constrain all vessel vertices?
        if (fixallvesselvertices == 1) then
            do i = 1, size(grid%bnd, 1)
                select case (grid%bnd(i)%ID)

                case (targetID, vesselID)

                    ! Get boundary vertices
                    tv = grid%bnd(i)%vert

                    ! Set mask
                    mask(tv) = .true.
                    do k = 1, nv
                        isnofieldline(k) = (grid%vert%fieldlineID(k) == 0_I8)
                    end do
                    do k = 1, nv
                        istempfree(k) = (cc(k) >= maxcc(k))
                    end do
                    ! De Morgan: .not.(A .or. istp .or. isconstrainedv)
                    !          = (.not. A) .and. (.not. istp) .and. (.not. isconstrainedv)
                    mask = mask .and. (.not. istempfree) .and. (.not. istp) &
                        .and. (.not. isconstrainedv) .and. isnofieldline
                    where (mask) isfixedpoint = .true.
                    deallocate(tv)

                case default

                    ! Do nothing

                end select

            end do
        end if

        ! Constrain the endpoints of target plates to their own flux
        ! values, in order to avoid shitty behaviour when having
        ! targets that are nearly flux aligned.
        if (fixfluxalignedtargets == 1) then
            do i = 1, size(grid%bnd, 1)
                select case (grid%bnd(i)%ID)

                case (targetID)

                    ! Get the end vertices
                    tv = grid%Bnd(i)%vert([1, size(grid%Bnd(i)%vert, 1)])

                    ! Delete vertices that can't be constrained or are
                    ! constrained already
                    mask(tv) = .true.
                    do k = 1, nv
                        istempfree(k) = (cc(k) >= maxcc(k))
                    end do
                    where (isconstrainedv .or. istempfree .or. istp) mask = .false.
                    where (mask) isfixedpoint = .true.
                    deallocate(tv)

                case default

                    ! Do nothing

                end select
            end do
        end if

        ! Constrain the points on the vessel outermost boundary, if it
        ! exists. The points are only constrained if there exists a
        ! non-boundary vertex that is a neighbour of the current
        ! vertex and which has the same flux line ID as the current
        ! vertex, which should be nonzero. Otherwise, the flux value is
        ! NOT constrained.
        if (fixfarvesselflux == 1) then
            do i = 1, size(grid%bnd, 1)
                select case (grid%bnd(i)%ID)

                case default

                    ! Do nothing

                case (vesselID)

                    ! Get the vertices
                    tv = grid%Bnd(i)%vert
                    where (vert%fieldlineID(tv) == 0) mask(tv) = .true.
                    deallocate(tv)

                    ! Delete vertices that can't be constrained or are
                    ! constrained already
                    do k = 1, nv
                        istempfree(k) = (cc(k) >= maxcc(k))
                    end do
                    where ( istempfree .or. (isconstrainedv) .or. (istp)) mask = .false.

                end select
            end do
        end if

        ! Override if desired - include vertices that have not been
        ! constrained as fixed points
        if (options%ffoptions%doboxoverride == 1) then
            ! Determine which points are already constrained - these
            ! are all points occuring in fixed points, vpairs, and
            ! x-points.
            allocate(includev(nv))
            includev(:) = .false.
            boxx = options%ffoptions%includeboxx
            boxy = options%ffoptions%includeboxy
            if (size(boxx, 1) > 0) then
                ! Loop over all boxes to include
                do k = 1, size(boxx, 1)
                    where ((vert%x > boxx(k, 1)) .and. (vert%x < boxx(k, 2)) &
                        .and. (vert%y > boxy(k, 1)) .and. (vert%y < boxy(k, 2))) &
                        includev = .true.
                end do
            end if

            ! Add points that have not yet been constrained
            do k = 1, nv
                istempfree(k) = (cc(k) >= maxcc(k))
            end do
            where (isconstrainedv .or. istempfree) includev = .false.
            where (includev) mask = .true.

        end if

        ! Add fixed points
        if (count(mask) > 0) then
            allocate(tv(count(mask)))
            tv = pack(vID, mask)
            nfp = size(tv, 1)
            allocate(constraints%fixedpoints(nfp))
            do i = 1, nfp
                ! Add as fixed point
                allocate(constraints%fixedpoints(i)%ID(1))
                constraints%fixedpoints(i)%ID(1) = tv(i)
                constraints%fixedpoints(i)%PsiD = PsiD_tmp(tv(i))
                constraints%fixedpoints(i)%nID = 1
                constraints%fixedpoints(i)%fsID = fieldlineID(tv(i))
            end do

            ! Update counter
            isconstrainedv(tv) = .true.
            cc(tv) = cc(tv) + 1

            ! Deallocate
            deallocate(tv)
        end if

        ! Classical flux surfaces
        !========================
        ! Initialize
        nfs = size(fsIDs, 1)
        allocate(constraints%fluxsurfaces(nfs))

        ! Determine vessel vertices (for initial PsiD later on)
        call DetermineVesselVertices(isvesselvertex, isvesselface, grid)

        ! Loop over all the flux surfaces to compute desired flux
        do i = 1, nfs

            ! Get all vertices with this ID
            mask(:) = (grid%vert%fieldlineID == fsIDs(i))

            ! Get the flux values
            if (any(pack(isvesselvertex, mask))) then
                ! Average only over boundary vertices
                tpsi = sum(pack(PsiD_tmp, (mask .and. isvesselvertex))) &
                    /count((mask .and. isvesselvertex))
            else
                ! Average over all vertices
                tpsi = sum(pack(PsiD_tmp, mask))/count(mask)
            end if

            ! Add value to flux surface
            constraints%fluxsurfaces(i)%PsiD    = tpsi

            ! Check which IDs can be added and add them
            do k = 1, nv
                istempfree(k) = (cc(k) >= maxcc(k))
            end do
            mask = mask .and. (  ( (.not. isconstrainedv) .or. &
                (.not. istempfree) ) .or. ( isxp .or. istp ) ) ! keep x-points and tps

            ! Add
            allocate(constraints%fluxsurfaces(i)%ID(count(mask)))
            constraints%fluxsurfaces(i)%ID      = pack(vID, mask)
            constraints%fluxsurfaces(i)%nID     = size(constraints%fluxsurfaces(i)%ID, 1)
            constraints%fluxsurfaces(i)%PsiD    = tpsi
            constraints%fluxsurfaces(i)%fsID    = fsIDs(i)

            ! Update counter
            cc(constraints%fluxsurfaces(i)%ID) = cc(constraints%fluxsurfaces(i)%ID) + 1
            isconstrainedv(constraints%fluxsurfaces(i)%ID) = .true.

        end do

        ! Compute total number of constraints
        !====================================
        constraints%ncon = nfp + ntp
        nfscon = 0
        nspcon = 0
        nfpcon = nfp
        ntpcon = ntp
        do i = 1, nfs
            constraints%ncon = constraints%ncon + constraints%fluxsurfaces(i)%nID
            nfscon = nfscon + constraints%fluxsurfaces(i)%nID
        end do
        do i = 1, nsp
            constraints%ncon = constraints%ncon + constraints%specialpoints(i)%nID-1
            nspcon = nspcon + constraints%specialpoints(i)%nID-1
        end do

        ! Debugging info
        !===============
        ! Write datafile
        if (options%writedata == 1) then
            call constraints%WriteData(grid)
        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateFluxfunctionConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! Evaluate the flux function constraints imposed on the
        ! vertices. For each  vertex considered (see InitDesign), the
        ! flux function is imposed mathematically as:
        !
        !       G_i = Psi(x_i,y_i) - Psi_D,
        !
        ! where x_i and y_i are the i-th vertex's coordinates, Psi is
        ! the underlying flux function, characterized by a bicubic
        ! spline interpolant, and Psi_D is a vector containing the
        ! desired flux function values.

        ! The Hessian of each ith constraint is:
        !
        !       Hjk,i = d^2F/dx^2 if j == k == i
        !       Hjk,i = d^2F/dy^2 if j == k == i+numel(x)
        !       Hjk,i = d^2F/dxdy if j == k+numel(x)
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = d^2F/dx^2 lambda_i (idem for y)

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is
        ! returned, but the hessian-vector multiplication with the
        ! vector lambda, which should be of suitable size.

        ! Note 2: the row and column indices for the constraints are
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine.

        ! Note 3: at first, we compute the linearization of the
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e.
        ! transpose) to obtain the gradient.

        ! Initialize
        !===========
        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(FluxfunctionConstraintsUDT)   :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, &
            hessG_coord, gradG_coord, hessG_flux, gradG_flux
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, i

        ! Auxiliary variables
        real(R8), allocatable               :: psival(:), dpsidx(:), &
            dpsidy(:), d2psidx2(:), d2psidxdy(:), d2psidy2(:), &
            dVdx(:), dVdy(:), d2Vdx2(:), d2Vdxdy(:), d2Vdy2(:), &
            d3Vdx3(:), d3Vdx2dy(:), d3Vdxdy2(:), d3Vdy3(:), &
            d3psidx3(:), d3psidx2dy(:), d3psidxdy2(:), d3psidy3(:), &
            Gv(:)
        integer(I8)                         :: spID, nc
        integer(I8), allocatable            :: tvID(:)

        ! Data

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if
        if (present(dGdvarin)) then
            dGdvar = dGdvarin
        end if
        if (present(dgradGdvarin)) then
            dgradGdvar = dgradGdvarin
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if
        if (allocated(G)) then
            if (size(G, 1) .ne. constraints%ncon) then
                deallocate(G)
                allocate(G(constraints%ncon))
            end if
        else
            allocate(G(constraints%ncon))
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            plf     => environment%vessel%plfvessel,                &
            nv      => grid%vert%ntot,                              &
            specialpoints           => constraints%specialpoints,   &
            nsp                     => constraints%nspecialpoints,  &
            fluxsurfaces            => constraints%fluxsurfaces,    &
            nfs                     => constraints%nfluxsurfaces,   &
            fixedpoints             => constraints%fixedpoints,     &
            nfp                     => constraints%nfixedpoints,    &
            tangencypoints          => constraints%tangencypoints,  &
            ntp                     => constraints%ntangencypoints, &
            interp  => magneticField%interp,    &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Pre-evaluate some data
        allocate(psival(nv))
        allocate(d3psidx3(nv), d3psidx2dy(nv), d3psidxdy2(nv), &
            d3psidy3(nv), d3Vdx3(nv), d3Vdx2dy(nv), d3Vdxdy2(nv), d3Vdy3(nv))
        allocate(d2Vdx2(nv), d2Vdxdy(nv), d2Vdy2(nv))
        allocate(d2psidx2(nv), d2psidxdy(nv), d2psidy2(nv))
        allocate(dpsidx(nv), dpsidy(nv))
        allocate(dVdx(nv), dVdy(nv))
        allocate(Gv(nv))
        Gv = 0.0_R8

        call magneticField%interp%Evaluate(x, y, 0, 0, psival)

        if (ntp > 0) then
            ! Evaluate vessel shape derivatives
            call plf%Evaluate(x, y, 1, 0, dVdx)
            call plf%Evaluate(x, y, 0, 1, dVdy)

            if (.not. dogradient) then
                ! Precompute
                call magneticField%interp%Evaluate(x, y, 1, 0, dpsidx)
                call magneticField%interp%Evaluate(x, y, 0, 1, dpsidy)
            end if
            if ( (.not. dohessian) .and. dogradient) then
                ! Precompute
                call magneticField%interp%Evaluate(x, y, 2, 0, d2psidx2)
                call magneticField%interp%Evaluate(x, y, 1, 1, d2psidxdy)
                call magneticField%interp%Evaluate(x, y, 0, 2, d2psidy2)
            end if

            ! Evaluate additional derivatives
            if (dogradient) then
                call plf%Evaluate(x, y, 2, 0, d2Vdx2)
                call plf%Evaluate(x, y, 1, 1, d2Vdxdy)
                call plf%Evaluate(x, y, 0, 2, d2Vdy2)
            end if
            if (dohessian) then
                ! Precompute
                call magneticField%interp%Evaluate(x, y, 3, 0, d3psidx3)
                call magneticField%interp%Evaluate(x, y, 2, 1, d3psidx2dy)
                call magneticField%interp%Evaluate(x, y, 1, 2, d3psidxdy2)
                call magneticField%interp%Evaluate(x, y, 0, 3, d3psidy3)

                call plf%Evaluate(x, y, 3, 0, d3Vdx3)
                call plf%Evaluate(x, y, 2, 1, d3Vdx2dy)
                call plf%Evaluate(x, y, 1, 2, d3Vdxdy2)
                call plf%Evaluate(x, y, 0, 3, d3Vdy3)
            end if
        end if


        if (dogradient) then
            ! Precompute
            call magneticField%interp%Evaluate(x, y, 1, 0, dpsidx)
            call magneticField%interp%Evaluate(x, y, 0, 1, dpsidy)
        end if

        if (dohessian) then
            ! Precompute
            call magneticField%interp%Evaluate(x, y, 2, 0, d2psidx2)
            call magneticField%interp%Evaluate(x, y, 1, 1, d2psidxdy)
            call magneticField%interp%Evaluate(x, y, 0, 2, d2psidy2)
        end if

        !if (ntp > 0) then
        !    ! Compute vessel boundary derivatives etc as well
        !    ! For now, call error
        !    call gdErrorHandler('Tangency points not yet implemented')
        !end if

        ! Evaluate
        !=========
        ic = 0
        ! Special points
        do i = 1, nsp
            ! Unpack
            nc = specialpoints(i)%nID-1
            spID = specialpoints(i)%ID(1)
            tvID = specialpoints(i)%ID(2:specialpoints(i)%nID)

            ! Evaluate
            G(ic+1:ic+nc) = psival(tvID) - psival(spID)
            Gv(tvID) = G(ic+1:ic+nc)

            ! Update
            ic = ic + nc
        end do

        ! Tangency points
        do i = 1, ntp
            ! Unpack
            nc = tangencypoints(i)%nID
            tvID = tangencypoints(i)%ID

            ! Evaluate
            G(ic+1:ic+nc) = -dpsidx(tvID)*dVdy(tvID) + dpsidy(tvID)*dVdx(tvID)
            Gv(tvID) = G(ic+1:ic+nc)

            ! Update
            ic = ic + nc
        end do

        ! Fixed points
        do i = 1, nfp
            ! Unpack
            nc = fixedpoints(i)%nID
            tvID = fixedpoints(i)%ID

            ! Evaluate
            G(ic+1:ic+nc) = psival(tvID) - fixedpoints(i)%PsiD
            Gv(tvID) = G(ic+1:ic+nc)

            ! Update
            ic = ic + nc
        end do

        ! Flux surfaces
        do i = 1, nfs
            ! Unpack
            nc = fluxsurfaces(i)%nID
            tvID = fluxsurfaces(i)%ID

            ! Evaluate
            G(ic+1:ic+nc) = psival(tvID) - fluxsurfaces(i)%PsiD
            Gv(tvID) = G(ic+1:ic+nc)

            ! Update
            ic = ic + nc
        end do

        ! write to file
        !call Write3DCoordinateData(grid%vert%x, grid%vert%y, Gv, 'con_ff_val_vertices')

        ! Derivatives
        !============
        ! Set row and column sizes
        gradG%nrow = designvariables%nphi
        gradG%ncol = constraints%ncon
        hessG%nrow = designvariables%nphi
        hessG%ncol = designvariables%nphi

        ! Check which gradient to compute
        select case (designvariables%type)

        case ('coordinates')

            call constraints%EvaluateDerivativesCoordinates(grid, gradG, &
                hessG, dogradient, dohessian, lambda, psival, dpsidx, &
                dpsidy, d2psidx2, d2psidxdy, d2psidy2, d3psidx3, &
                d3psidx2dy, d3psidxdy2, d3psidy3, dVdx, dVdy, d2Vdx2, &
                d2Vdxdy, d2Vdy2, d3Vdx3, d3Vdx2dy, d3Vdxdy2, d3Vdy3)

        case ('coordinates_desiredflux')

            ! Coordinate contribution
            gradG_coord%nrow = designvariables%nphi
            gradG_coord%ncol = constraints%ncon
            hessG_coord%nrow = designvariables%nphi
            hessG_coord%ncol = designvariables%nphi
            call constraints%EvaluateDerivativesCoordinates(grid, gradG_coord, &
                hessG_coord, dogradient, dohessian, lambda, psival, dpsidx, &
                dpsidy, d2psidx2, d2psidxdy, d2psidy2, d3psidx3, &
                d3psidx2dy, d3psidxdy2, d3psidy3, dVdx, dVdy, d2Vdx2, &
                d2Vdxdy, d2Vdy2, d3Vdx3, d3Vdx2dy, d3Vdxdy2, d3Vdy3)

            ! Flux contribution
            gradG_flux%nrow = designvariables%nphi
            gradG_flux%ncol = constraints%ncon
            hessG_flux%nrow = designvariables%nphi
            hessG_flux%ncol = designvariables%nphi
            call constraints%EvaluateDerivativesFlux(gradG_flux, hessG_flux, &
                dogradient, dohessian, lambda)

            ! Update flux contribution design variable indices
            gradG_flux%row = gradG_flux%row + 2*grid%vert%ntot
            hessG_flux%row = hessG_flux%row + 2*grid%vert%ntot
            hessG_flux%col = hessG_flux%col + 2*grid%vert%ntot

            ! Combine
            gradG = gradG_coord + gradG_flux
            hessG = hessG_coord + hessG_flux

        case default

            call gdErrorHandler('EvaluateFluxfunctionConstraints: unknown design variable type')

        end select

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    ! Gradient & Hessian computation, coordinates
    subroutine EvaluateCoordinatesDerivativesFluxFunctionConstraints(&
        constraints, grid, gradG, hessG, dogradient, dohessian, &
        lambda, psival, dpsidx, dpsidy, d2psidx2, d2psidxdy, d2psidy2, &
        d3psidx3, d3psidx2dy, d3psidxdy2, d3psidy3, dVdx, dVdy, &
        d2Vdx2, d2Vdxdy, d2Vdy2, d3Vdx3, d3Vdx2dy, d3Vdxdy2, d3Vdy3)

        ! Description
        !============
        ! This routine evaluates the gradient and hessian w.r.t. the
        ! grid coordinates. It is assumed that the number of rows and
        ! columns is already computed before.

        ! Note: the derivatives are computed as if the coordinates are
        ! the only design variables. Any reordering/adjustment of column
        ! or row indices should be done afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(FluxfunctionConstraintsUDT)   :: constraints
        real(R8), dimension(:), intent(in)  :: lambda(*), psival(*), &
            dpsidx(*), dpsidy(*), d2psidx2(*), d2psidxdy(*), d2psidy2(*), &
            d3psidx3(*), d3psidx2dy(*), d3psidxdy2(*), d3psidy3(*), &
            dVdx(*), dVdy(*), d2Vdx2(*), d2Vdxdy(*), d2Vdy2(*), &
            d3Vdx3(*), d3Vdx2dy(*), d3Vdxdy2(*), d3Vdy3(*)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT), intent(in)           :: grid
        logical                             :: dogradient, dohessian

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k, i
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        integer(I8)                         :: nvg, nvh, spID, nc
        integer(I8), allocatable            :: tvID(:)

        ! Initialize
        !===========
        ! Associate
        associate(&
            nv      => grid%vert%ntot,          &
            specialpoints           => constraints%specialpoints,   &
            nsp                     => constraints%nspecialpoints,  &
            fluxsurfaces            => constraints%fluxsurfaces,    &
            nfs                     => constraints%nfluxsurfaces,   &
            fixedpoints             => constraints%fixedpoints,     &
            nfp                     => constraints%nfixedpoints,    &
            tangencypoints          => constraints%tangencypoints,  &
            ntp                     => constraints%ntangencypoints, &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Counters
        ic = 0
        ivg = 0
        ivh = 0

        ! Check allocation
        ! Precompute number of entries for jacobians/hessian
        if (dogradient .or. dohessian) then
            ! Compute number of jacobian and hessian values
            nvg = 0
            nvh = 0

            ! Special point contributions
            do i = 1, nsp
                ! Gradient: 4 entries per non-special point
                nvg = nvg + 4*(specialpoints(i)%nID-1)

                ! Hessian: 8 entries per non-special point (but
                ! distributed over valxx, valxy etc)
                nvh = nvh + 8*(specialpoints(i)%nID-1)
            end do

            ! Fixed point contributions
            do i = 1, nfp
                ! Gradient: 2 entries per fixed point
                nvg = nvg + 2*fixedpoints(i)%nID

                ! Hessian: 4 entries per fixed point
                nvh = nvh + 4*fixedpoints(i)%nID
            end do

            ! Tangency point contributions
            do i = 1, ntp
                ! Gradient: 2 entries per tangency point
                nvg = nvg + 2*tangencypoints(i)%nID

                ! Hessian: 4 entries per tangency point
                nvh = nvh + 4*tangencypoints(i)%nID
            end do

            ! Flux surface contributions
            do i = 1, nfs
                ! Gradient: 2 entries per contribution
                nvg = nvg + 2*fluxsurfaces(i)%nID

                ! Hessian: 4 entires per contribution
                nvh = nvh + 4*fluxsurfaces(i)%nID
            end do

            if (.not. allocated(jacG%row)) then
                jacG%nval = nvg
                call jacG%Allocate()
            end if
            if (.not. allocated(hessG%row)) then
                hessG%nval = nvh
                call hessG%Allocate()
            end if
        else
            ! Nothing to compute
            if (.not. allocated(jacG%row)) then
                jacG%nval = 0
                call jacG%Allocate()
            end if
            if (.not. allocated(hessG%row)) then
                hessG%nval = 0
                call hessG%Allocate()
            end if
        end if

        ! Special points
        !===============
        do i = 1, nsp
            ! Unpack
            nc = specialpoints(i)%nID-1
            spID = specialpoints(i)%ID(1)
            tvID = specialpoints(i)%ID(2:specialpoints(i)%nID)
            conindex = [(k, k = ic+1, ic+nc)]

            ! Gradient
            if (dogradient) then
                ! x
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID
                jacG%val(valindex) = dpsidx(tvID)
                ivg = ivg + nc

                valindex = valindex + nc
                jacG%row(valindex) = conindex
                jacG%col(valindex) = spID
                jacG%val(valindex) = -dpsidx(spID)
                ivg = ivg + nc

                ! y
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID + grid%vert%ntot
                jacG%val(valindex) = dpsidy(tvID)
                ivg = ivg + nc

                valindex = valindex + nc
                jacG%row(valindex) = conindex
                jacG%col(valindex) = spID + grid%vert%ntot
                jacG%val(valindex) = -dpsidy(spID)
                ivg = ivg + nc

            end if

            ! Hessian
            if (dohessian) then
                ! xx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidx2(tvID)*lambda(conindex)
                ivh = ivh + nc

                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = spID
                hessG%col(valindex) = spID
                hessG%val(valindex) = -d2psidx2(spID)*lambda(conindex)
                ivh = ivh + nc

                ! xy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = spID
                hessG%col(valindex) = spID + grid%vert%ntot
                hessG%val(valindex) = -d2psidxdy(spID)*lambda(conindex)
                ivh = ivh + nc

                ! yx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = spID + grid%vert%ntot
                hessG%col(valindex) = spID
                hessG%val(valindex) = -d2psidxdy(spID)*lambda(conindex)
                ivh = ivh + nc

                ! yy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidy2(tvID)*lambda(conindex)
                ivh = ivh + nc

                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = spID + grid%vert%ntot
                hessG%col(valindex) = spID + grid%vert%ntot
                hessG%val(valindex) = -d2psidy2(spID)*lambda(conindex)
                ivh = ivh + nc

            end if

            ! Update counter
            ic = ic + nc
        end do

        ! Tangency points
        !================
        ! Not yet implemented
        do i = 1, ntp
            ! Unpack
            nc = tangencypoints(i)%nID
            tvID = tangencypoints(i)%ID
            conindex = [(k, k = ic+1, ic+nc)]

            ! Gradient
            if (dogradient) then
                ! x
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID
                jacG%val(valindex) = -d2psidx2(tvID)*dVdy(tvID) - &
                    dpsidx(tvID)*d2Vdxdy(tvID) + d2psidxdy(tvID)*dVdx(tvID) &
                    + dpsidy(tvID)*d2Vdx2(tvID)
                ivg = ivg + nc

                ! y
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID + grid%vert%ntot
                jacG%val(valindex) = -d2psidxdy(tvID)*dVdy(tvID) - &
                    dpsidx(tvID)*d2Vdy2(tvID) + d2psidy2(tvID)*dVdx(tvID) &
                    + dpsidy(tvID)*d2Vdxdy(tvID)
                ivg = ivg + nc
            end if

            ! Hessian
            if (dohessian) then
                ! xx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID
                hessG%val(valindex) = (-d3psidx3(tvID)*dVdy(tvID) &
                    - d2psidx2(tvID)*d2Vdxdy(tvID) - d2psidx2(tvID)*d2Vdxdy(tvID) &
                    - dpsidx(tvID)*d3Vdx2dy(tvID) + d3psidx2dy(tvID)*dVdx(tvID) &
                    + d2psidxdy(tvID)*d2Vdx2(tvID) + d2psidxdy(tvID)*d2Vdx2(tvID) &
                    + dpsidy(tvID)*d3Vdx3(tvID))*lambda(conindex)
                ivh = ivh + nc

                ! xy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = (-d3psidx2dy(tvID)*dVdy(tvID) &
                    - d2psidx2(tvID)*d2Vdy2(tvID) - d2psidxdy(tvID)*d2Vdxdy(tvID) &
                    - dpsidx(tvID)*d3Vdxdy2(tvID) + d3psidxdy2(tvID)*dVdx(tvID) &
                    + d2psidxdy(tvID)*d2Vdxdy(tvID) + d2psidy2(tvID)*d2Vdx2(tvID) &
                    + dpsidy(tvID)*d3Vdx2dy(tvID))*lambda(conindex)
                ivh = ivh + nc

                ! yx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID
                hessG%val(valindex) = (-d3psidx2dy(tvID)*dVdy(tvID) &
                    - d2psidx2(tvID)*d2Vdy2(tvID) - d2psidxdy(tvID)*d2Vdxdy(tvID) &
                    - dpsidx(tvID)*d3Vdxdy2(tvID) + d3psidxdy2(tvID)*dVdx(tvID) &
                    + d2psidxdy(tvID)*d2Vdxdy(tvID) + d2psidy2(tvID)*d2Vdx2(tvID) &
                    + dpsidy(tvID)*d3Vdx2dy(tvID))*lambda(conindex)
                ivh = ivh + nc

                ! yy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = (-d3psidxdy2(tvID)*dVdy(tvID) &
                    - d2psidxdy(tvID)*d2Vdy2(tvID) - d2psidxdy(tvID)*d2Vdy2(tvID) &
                    - dpsidx(tvID)*d3Vdy3(tvID) + d3psidy3(tvID)*dVdx(tvID) &
                    + d2psidy2(tvID)*d2Vdxdy(tvID) + d2psidy2(tvID)*d2Vdxdy(tvID) &
                    + dpsidy(tvID)*d3Vdxdy2(tvID))*lambda(conindex)
                ivh = ivh + nc

            end if

            ! Update counter
            ic = ic + nc
        end do

        ! Fixed points
        !=============
        do i = 1, nfp
            ! Unpack
            nc = fixedpoints(i)%nID
            tvID = fixedpoints(i)%ID
            conindex = [(k, k = ic+1, ic+nc)]

            ! Gradient
            if (dogradient) then
                ! x
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID
                jacG%val(valindex) = dpsidx(tvID)
                ivg = ivg + nc

                ! y
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID + grid%vert%ntot
                jacG%val(valindex) = dpsidy(tvID)
                ivg = ivg + nc
            end if

            ! Hessian
            if (dohessian) then
                ! xx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidx2(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! xy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! yx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! yy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidy2(tvID)*lambda(conindex)
                ivh = ivh + nc

            end if

            ! Update counter
            ic = ic + nc
        end do

        ! Flux surfaces
        !==============
        do i = 1, nfs
            ! Unpack
            nc = fluxsurfaces(i)%nID
            tvID = fluxsurfaces(i)%ID
            conindex = [(k, k = ic+1, ic+nc)]

            ! Gradient
            if (dogradient) then
                ! x
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID
                jacG%val(valindex) = dpsidx(tvID)
                ivg = ivg + nc

                ! y
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tvID + grid%vert%ntot
                jacG%val(valindex) = dpsidy(tvID)
                ivg = ivg + nc
            end if

            ! Hessian
            if (dohessian) then
                ! xx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidx2(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! xy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! yx
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID
                hessG%val(valindex) = d2psidxdy(tvID)*lambda(conindex)
                ivh = ivh + nc

                ! yy
                valindex = [(k, k = ivh+1, ivh+nc)]
                hessG%row(valindex) = tvID + grid%vert%ntot
                hessG%col(valindex) = tvID + grid%vert%ntot
                hessG%val(valindex) = d2psidy2(tvID)*lambda(conindex)
                ivh = ivh + nc

            end if

            ! Update counter
            ic = ic + nc

        end do


        ! Transpose
        !==========
        gradG%row = jacG%col
        gradG%col = jacG%row
        gradG%val = jacG%val
        gradG%nval = jacG%nval

        ! Housekeeping
        end associate


    end subroutine

    ! Gradient & Hessian computation, psi values
    subroutine EvaluateFluxDerivativesFluxFunctionConstraints(&
        constraints, gradG, hessG, dogradient, dohessian, &
        lambda)

        ! Description
        !============
        ! This routine evaluates the gradient and hessian w.r.t. the
        ! desired flux values. It is assumed that the number of rows and
        ! columns is already computed before.

        ! Note: the derivatives are computed as if the coordinates are
        ! the only design variables. Any reordering/adjustment of column
        ! or row indices should be done afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(FluxfunctionConstraintsUDT)       :: constraints
        type(MySparseUDT)                       :: gradG, hessG
        logical, intent(in)                     :: dogradient, dohessian
        real(R8),  intent(in)                   :: lambda(*)

        ! Auxiliary
        type(MySparseUDT)                       :: jacG
        integer(I8)                             :: nc

        ! Loop
        integer(I8)                             :: i, k, ic

        ! Compute
        !========
        ! Derivatives are simply equal to minus one, just need to
        ! correctly determine columns and rows.
        ! Jacobian
        jacG%nval = constraints%nfscon
        call jacG%Allocate()

        ! Columns are design variables for Jacobian, so simply equal
        ! to flux function constraint index (locally)
        ic = 0
        do i = 1, constraints%nfluxsurfaces
            ! Get number of constraints
            nc = constraints%fluxsurfaces(i)%nID

            ! Set
            jacG%col(ic+1:ic+nc) = i

            ! Update
            ic = ic + nc
        end do

        ! Rows are constraint indices - flux function constraints are
        ! evaluated last
        ic = constraints%nfpcon + constraints%ntpcon + constraints%nspcon
        jacG%row = [(k, k = ic+1, ic+constraints%nfscon)]

        ! Values are just one
        jacG%val(:) = -1

        ! Compute gradient
        !=================
        gradG%val = jacG%val
        gradG%row = jacG%col
        gradG%col = jacG%row
        gradG%nval = jacG%nval

        ! Hessian
        !========
        ! Simply zero
        hessG%nval = 0
        call hessG%Allocate()

    end subroutine

    ! Data output
    subroutine WriteDataFluxfunctionConstraints(constraints, grid)

        ! Description
        !============
        ! Write grid nodes in the following format:
        ! ID, x, y
        ! Different files are written for special vertices, fixed
        ! vertices, flux surface vertices, and tangency points

        ! The usual
        implicit none

        ! Declare variables
        type(gridUDT), intent(in)                       :: grid
        class(FluxfunctionConstraintsUDT)               :: constraints
        real(R8), allocatable                           :: x(:), y(:)
        integer(I8)                                     :: nIDs, nIDsfs, &
            counter, i
        integer(I8), allocatable                        :: IDs(:)
        character(:), allocatable                       :: filepath, &
            thispath

        ! Initialize
        !===========
        ! Set the correct directories
        allocate(character(len('con_ff_vertices')) :: filepath)
        filepath = 'con_ff_vertices'

        ! Associate
        associate(&
            nsp             => constraints%nspecialpoints,  &
            sp              => constraints%specialpoints,   &
            ntp             => constraints%ntangencypoints, &
            tp              => constraints%tangencypoints,  &
            nfp             => constraints%nfixedpoints,    &
            fp              => constraints%fixedpoints,     &
            nfs             => constraints%nfluxsurfaces,   &
            fs              => constraints%fluxsurfaces     &
            )

        ! Special points
        !===============
        ! Extract IDs of special points only
        nIDs    = nsp
        nIDsfs  = 0 ! %already compute how much special points we will have
        allocate(IDs(nIDs), x(nIDs), y(nIDs))
        do i = 1, nsp
            IDs(i) = sp(i)%ID(1)
            nIDsfs = nIDsfs + sp(i)%nID-1
        end do
        x = grid%vert%x(IDs)
        y = grid%vert%y(IDs)

        ! Write
        thispath = filepath // '_sp'
        call WriteVertexData(IDs, x, y, thispath)

        ! Deallocate
        deallocate(x ,y, IDs)

        ! Extract IDs of flux surface vertices with special points (but without those points)
        allocate(IDs(nIDsfs), x(nIDsfs),y(nIDsfs))
        counter = 0
        do i = 1, nsp
            ! Get number of IDs
            nIDs = sp(i)%nID

            ! Add coordinates
            IDs(counter+1:counter+nIDs-1) = sp(i)%ID(2:nIDs)
            x(counter+1:counter+nIDs-1) = grid%vert%x(sp(i)%ID(2:nIDs))
            y(counter+1:counter+nIDs-1) = grid%vert%y(sp(i)%ID(2:nIDs))

            ! Update counter
            counter = counter + nIDs-1
        end do

        ! Write
        thispath = filepath // '_spfs'
        call WriteVertexData(IDs, x, y, thispath)

        ! Deallocate
        deallocate(x ,y, IDs)

        ! Fixed points
        !=============
        ! Extract IDs
        nIDs    = nfp
        allocate(IDs(nIDs), x(nIDs), y(nIDs))
        do i = 1, nfp
            IDs(i) = fp(i)%ID(1)
        end do
        x = grid%vert%x(IDs)
        y = grid%vert%y(IDs)

        ! Write
        thispath = filepath // '_fp'
        call WriteVertexData(IDs, x, y, thispath)

        ! Deallocate
        deallocate(x, y, IDs)

        ! Tangency points
        !================
        ! Extract IDs
        nIDs    = ntp
        allocate(IDs(nIDs), x(nIDs), y(nIDs))
        do i = 1, ntp
            IDs(i) = tp(i)%ID(1)
        end do
        x = grid%vert%x(IDs)
        y = grid%vert%y(IDs)

        ! Write
        thispath = filepath // '_tp'
        call WriteVertexData(IDs, x, y, thispath)

        ! Deallocate
        deallocate(x ,y, IDs)

        ! Flux surfaces
        !==============
        ! Compute total number of points
        nIDs    = 0
        do i = 1, nfs
            nIDs = nIDs + fs(i)%nID
        end do
        allocate(IDs(nIDs), x(nIDs), y(nIDs))

        ! Set points
        counter = 0
        do i = 1, nfs
            ! Get number of IDs
            nIDs = fs(i)%nID

            ! Add coordinates
            IDs(counter+1:counter+nIDs) = fs(i)%ID
            x(counter+1:counter+nIDs) = grid%vert%x(fs(i)%ID)
            y(counter+1:counter+nIDs) = grid%vert%y(fs(i)%ID)

            ! Update counter
            counter = counter + nIDs
        end do


        ! Write
        thispath = filepath // '_fs'
        call WriteVertexData(IDs, x, y, thispath)

        ! Deallocate
        deallocate(x ,y, IDs)

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Destructor
    subroutine DestroyFluxfunctionConstraints(constraints)

        ! Description
        !============
        ! Destructor of the flux function constraints

        ! Declare variables
        !==================
        ! Arguments
        type(FluxfunctionConstraintsUDT)        :: constraints

        ! Destroy
        !========

    end subroutine



    !------------------------------------------------------------------!
    !                         BOUNDARY FUNCTION                        !
    !------------------------------------------------------------------!
    ! Initialize
    subroutine InitializeBoundaryFunctionConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the boundary function constraints. Here, the
        ! function equals zero when a point lies precisely on the
        ! boundary. Though the 'boundary' on which the nodes should lie
        ! can be defined in many different ways, we chose here to base
        ! it on the vessel bounding polygon, which basically
        ! encompasses the plasma chamber.

        ! Only the vertices lying on the target plate are considered
        ! here to be constrained, though one could, if for some reason
        ! desired, also extend this towards other boundaries/boundary
        ! descriptions and other nodes.

        ! Modules
        use mod_definitions, only : targetID, vesselID

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)       :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                 :: designvariables

        ! Auxiliary
        logical                                     :: debugplots
        integer(I8)                                 :: ic, nv
        integer(I8), allocatable                    :: tv(:)
        logical, allocatable                        :: mask(:), &
            isconstrained(:)

        ! Loop
        integer(I8)                                 :: i, j

        ! Data
        data debugplots /.false./

        ! Initialize
        !===========
        ! Check allocation
        if (allocated(constraints%vert)) then
            deallocate(constraints%vert)
        end if

        ! Associate
        associate(&
            vessel      => environment%vessel  &
            )

        ! Bookkeeping of constrained vertices (to prevent imposing
        ! constraint twice)
        allocate(isconstrained(grid%vert%ntot))
        isconstrained = .false.

        ! Construct boundary
        !===================
        ! Should already be constructed in vessel - assign
        constraints%plf = vessel%plfvessel

        ! Visualize
        call constraints%plf%Visualize('constraints_boundary_plf')

        ! Set the constraints
        !====================
        ! Compute number of constraints
        constraints%ncon = 0
        constraints%nvert = 0
        do i = 1, size(grid%bnd)

            ! Check if target plate - hard coded here...
            if (any(grid%bnd(i)%ID == [targetID, vesselID])) then

                ! Get the current vertices
                allocate(tv(grid%bnd(i)%nvert))
                tv(:) = grid%bnd(i)%vert

                ! Construct the mask
                allocate(mask(grid%bnd(i)%nvert))
                mask(:) = .true.

                ! Check the monitor
                where (monitor%eqvcc(tv) .ge. monitor%maxeqvcc(tv)) mask = .false.

                ! Check if already constrained (will happen for nodes
                ! belonging to multiple boundaries or to boundaries that
                ! are closed upon themselves)
                do j = 1, size(tv)
                    if (isconstrained(tv(j))) then
                        mask(j) = .false.
                    else
                        isconstrained(tv(j)) = .true.
                    end if
                end do
                nv = count(mask)

                ! Add these nodes
                constraints%ncon = constraints%ncon + nv
                constraints%nvert = constraints%nvert + nv

                ! Deallocate
                deallocate(tv, mask)

            end if

        end do

        ! Allocate
        allocate(constraints%vert(constraints%nvert))

        ! Add vertices
        isconstrained = .false.
        ic = 0
        do i = 1, size(grid%bnd)
            if (any(grid%bnd(i)%ID == [targetID, vesselID])) then

                ! Get the current vertices
                allocate(tv(grid%bnd(i)%nvert))
                tv(:) = grid%bnd(i)%vert

                ! Construct the mask
                allocate(mask(grid%bnd(i)%nvert))
                mask(:) = .true.

                ! Check the monitor
                where (monitor%eqvcc(tv) .ge. monitor%maxeqvcc(tv)) mask = .false.

                ! Check if already constrained (will happen for nodes
                ! belonging to multiple boundaries or to boundaries that
                ! are closed upon themselves)
                do j = 1, size(tv)
                    if (isconstrained(tv(j))) then
                        mask(j) = .false.
                    else
                        isconstrained(tv(j)) = .true.
                    end if
                end do
                nv = count(mask)

                ! Add these nodes
                constraints%vert(ic+1:ic+nv) = pack(tv, mask)

                ! Update counter
                ic = ic + nv

                ! Update monitor
                monitor%eqvcc(pack(tv, mask)) = &
                    monitor%eqvcc(pack(tv, mask)) + 1

                ! Deallocate
                deallocate(tv, mask)

            end if

        end do

        ! End associate
        end associate

        ! Write datafile
        if (options%writedata == 1) then
            call WriteBoundaryConstraintVertices(grid, constraints%vert)
        end if

    end subroutine

    ! Evaluation
    subroutine EvaluateBoundaryFunctionConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! Evaluate the boundary function constraints based on the
        ! boundary function F. For each point, this is given as:
        !
        !       G_i = F(x_i,y_i) = 0
        !
        ! where x_i and y_i are the i-th vertex's coordinates, F is
        ! the underlying boundary shape function which is only zero
        ! on the boundary itself. F is given by a boundary shape
        ! function type in the constraints (currently it has to be
        ! the same for each node, though this may be extended in
        ! the future by extending the polygon structure and vertex
        ! IDs)

        ! The Hessian of each ith constraint is:
        !
        !       Hjk,i = d^2F/dx^2 if j == k == i
        !       Hjk,i = d^2F/dy^2 if j == k == i+numel(x)
        !       Hjk,i = d^2F/dxdy if j == k+numel(x)
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = d^2F/dx^2 lambda_i (idem for y)

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is
        ! returned, but the hessian-vector multiplication with the
        ! vector lambda, which should be of suitable size.

        ! Note 2: the row and column indices for the constraints are
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine.

        ! Note 3: at first, we compute the linearization of the
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e.
        ! transpose) to obtain the gradient.

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)   :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh

        ! Auxiliary variables
        real(R8), allocatable, dimension(:) :: Gv

        ! Initialize
        !===========
        ! Initialize vertex values
        allocate(Gv(grid%vert%ntot))
        Gv = 0.0_R8

        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            plf     => constraints%plf,         &
            tv      => constraints%vert,        &
            ntv     => constraints%nvert,       &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateBoundaryFunctionConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Evaluate
        select case (var)

        case ('no')

            ! No derivatives needed
            call plf%Evaluate(x(tv), y(tv), 0, 0, G)

        case ('vesselcoordinates')

            ! Derivatives needed w.r.t. polygon structure coordinates
            call plf%Evaluate(x(tv), y(tv), 0, 0, G, &
                'polygonsetcoordinates', values, dGdvar)

        end select
        !Gv(tv) = G

        ! Write data
        !call Write3DCoordinateData(grid%vert%x, grid%vert%y, Gv, 'con_bnd_val_vertices')

        ! Derivatives
        !============
        ! Set row and column sizes
        gradG%nrow = designvariables%nphi
        gradG%ncol = constraints%ncon
        hessG%nrow = designvariables%nphi
        hessG%ncol = designvariables%nphi

        ! Check which gradient to compute
        select case (designvariables%type)

        case ('coordinates', 'coordinates_desiredflux')

            ! Call dedicated routine - no flux contributions
            call constraints%EvaluateDerivativesCoordinates(gradG, &
                hessG, grid, dogradient, dohessian, lambda, designvariables, &
                var, values, dgradGdvar)

        case default

            call gdErrorHandler('EvaluateBoundaryFunctionConstraints: unknown design variable type')

        end select

        ! Housekeeping
        !=============
        end associate

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    ! Updateing
    subroutine UpdateBoundaryFunctionConstraints(constraints, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Update the boundary function description according to the
        ! given grid, magnetic field and environment. Can be used to
        ! update constraint parameters after external updating of these
        ! quantities (e.g. when doing shape optimization, vessel will
        ! change etc).

        ! Note: only boundary description here is updated.

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)       :: constraints
        type(GridUDT), intent(in)                   :: grid
        type(MagneticFieldUDT), intent(in)          :: magneticField
        type(EnvironmentUDT), intent(in)            :: environment

        ! Construct boundary
        !===================
        ! Should already be constructed in vessel - assign
        constraints%plf = environment%vessel%plfvessel

        ! Visualize
        !call constraints%plf%Visualize('constraints_boundary_plf')

    end subroutine

    ! Derivatives, coordinates
    subroutine EvaluateCoordinatesDerivativesBoundaryFunctionConstraints(&
        constraints, gradG, hessG, grid, dogradient, dohessian, lambda, &
        designvariables, var, values, dgradGdvar)

        ! Description
        !============
        ! This routine evaluates the gradient and hessian w.r.t. the
        ! grid coordinates. It is assumed that the number of rows and
        ! columns is already computed before.

        ! Note: the derivatives are computed as if the coordinates are
        ! the only design variables. Any reordering/adjustment of column
        ! or row indices should be done afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)   :: constraints
        class(DesignvariablesGDUDT)         :: designvariables
        real(R8), allocatable, intent(in)   :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT), intent(in)           :: grid
        logical, intent(in)                 :: dogradient, dohessian

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)
        integer(I8)                         :: nverttot
        integer(I8), allocatable            :: tvshift(:)

        ! Auxiliary variables
        real(R8), allocatable               :: dpsfdx(:), dpsfdy(:), &
            valxx(:), valxy(:), valyy(:)

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dgradGdvar, &
            dpsfdxdvar, dpsfdydvar, tempderivx, tempderivy

        ! Initialize
        !===========
        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            plf     => constraints%plf,         &
            tv      => constraints%vert,        &
            ntv     => constraints%nvert,       &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian
        nverttot = grid%vert%ntot
        allocate(tvshift(ntv))
        do k = 1, ntv
            tvshift(k) = tv(k) + nverttot
        end do

        ! Derivatives
        !============
        ! Gradient
        if (dogradient) then
            ! Allocate
            jacG%nval = 2*ntv
            call jacG%Allocate()
            allocate(dpsfdx(ntv))
            allocate(dpsfdy(ntv))
            allocate(conindex(ntv))
            allocate(valindex(ntv))

            ! Compute the derivative values
            select case (var)

            case ('no')

                ! No derivatives
                call plf%Evaluate(x(tv), y(tv), 1, 0, dpsfdx)
                call plf%Evaluate(x(tv), y(tv), 0, 1, dpsfdy)
                dpsfdxdvar = SpZeros(size(tv, 1), size(values, 1))
                dpsfdydvar = SpZeros(size(tv, 1), size(values, 1))

            case ('vesselcoordinates')

                ! Derivative w.r.t. polygonset vertices
                call plf%Evaluate(x(tv), y(tv), 1, 0, dpsfdx, &
                    'polygonsetcoordinates', values, dpsfdxdvar)
                call plf%Evaluate(x(tv), y(tv), 0, 1, dpsfdy, &
                    'polygonsetcoordinates', values, dpsfdydvar)

            case default

                ! Not implemented
                call gdErrorHandler('EvaluateCoordinatesDerivativesBoundaryFunctionConstraints: ' // &
                    'variable not implemented')

            end select

            ! x-contribution
            !---------------
            ! Build indices
            conindex = [(k, k = ic+1, ic+ntv)]
            valindex = [(k, k = ivg+1, ivg+ntv)]

            ! Add values
            jacG%row(valindex) = conindex
            jacG%col(valindex) = tv
            jacG%val(valindex) = dpsfdx

            ! y-contribution
            !---------------
            ! Build indices
            ivg = ivg + ntv
            valindex = valindex + ntv

            ! Add values
            jacG%row(valindex) = conindex
            jacG%col(valindex) = tvshift
            jacG%val(valindex) = dpsfdy

            ! Build gradient
            gradG%nval = jacG%nval

            call gradG%Allocate()
            gradG%row = jacG%col
            gradG%col = jacG%row
            gradG%val = jacG%val

            ! Derivatives
            !------------
            ! Full product between lambda and linearization
            tempderivx = lambda*dpsfdxdvar ! only local derivative, since rows go 1:ntv
            tempderivy = lambda*dpsfdydvar

            ! Expand
            tempderivx%nrow = designvariables%nphi
            tempderivx%ncol = size(values)
            tempderivx%row = tv(tempderivx%row)
            tempderivy%nrow = designvariables%nphi
            tempderivy%ncol = size(values)
            tempderivy%row = tv(tempderivy%row) + nverttot

            ! Add
            dgradGdvar = tempderivx + tempderivy

            ! Housekeeping
            call jacG%Deallocate()
        end if

        ! Hessian
        if (dohessian) then

            ! Allocate
            hessG%nval = 4*ntv
            if (.not. allocated(valindex)) then
                allocate(valindex(ntv))
            end if
            if (.not. allocated(conindex)) then
                allocate(conindex(ntv))
            end if
            if (.not. allocated(hessG%val)) then
                call hessG%Allocate()
            end if
            allocate(valxx(ntv))
            allocate(valxy(ntv))
            allocate(valyy(ntv))

            ! Initialize
            valxx(:) = 0
            valyy(:) = 0
            valxy(:) = 0

            ! Compute contributions
            call plf%Evaluate(x(tv), y(tv), 2, 0, valxx)
            call plf%Evaluate(x(tv), y(tv), 0, 2, valyy)
            call plf%Evaluate(x(tv), y(tv), 1, 1, valxy)

            ! xx-contribution
            !----------------
            k = 1
            ! Build indices
            valindex = [(k, k = ivh+1, ivh+ntv)]

            ! Add values
            hessG%row(valindex) = tv
            hessG%col(valindex) = tv
            hessG%val(valindex) = valxx*lambda ! element-wise mult.

            ! xy-contribution
            !----------------
            ! Build indices
            ivh = ivh + ntv
            valindex = valindex + ntv

            ! Add values
            hessG%row(valindex) = tv
            hessG%col(valindex) = tvshift
            hessG%val(valindex) = valxy*lambda ! element-wise mult.

            ! yx-contribution
            !----------------
            ! symmetric with xy
            ! Build indices
            ivh = ivh + ntv
            valindex = valindex + ntv

            ! Add values
            hessG%row(valindex) = tvshift
            hessG%col(valindex) = tv
            hessG%val(valindex) = valxy*lambda ! element-wise mult.

            ! yy-contribution
            !----------------
            ! Build indices
            ivh = ivh + ntv
            valindex = valindex + ntv

            ! Add values
            hessG%row(valindex) = tvshift
            hessG%col(valindex) = tvshift
            hessG%val(valindex) = valyy*lambda ! element-wise mult.

        end if

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Derivatives, flux
    subroutine EvaluateFluxDerivativesBoundaryFunctionConstraints(&
        constraints, gradG, hessG)

        ! Description
        !============
        ! This routine evaluates the gradient and hessian w.r.t. the
        ! grid coordinates. It is assumed that the number of rows and
        ! columns is already computed before.

        ! Note: the derivatives are computed as if the coordinates are
        ! the only design variables. Any reordering/adjustment of column
        ! or row indices should be done afterwards.

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)       :: constraints
        type(MySparseUDT)                           :: gradG, hessG

        ! Set values
        !===========
        ! Simply zero, no dependencies
        gradG%nval = 0
        call gradG%Allocate()

        hessG%nval = 0
        call hessG%Allocate()

    end subroutine

    ! Derivatives, vessel coordinates
    subroutine EvaluateVesselcoordinatesDerivativesBoundaryFunctionConstraints(&
        constraints, values, jacG, hessG, grid, magneticField, &
        environment, designvariables, lambda)

        ! Description
        !============
        ! Evaluate the jacobian (not gradient) of the constraints with
        ! respect to the vessel coordinates. Also add the linearization
        ! of the 'lambda*gradG' term as the 'hessian'. It is assumed that
        ! all fields are up to date.

        ! Declare variables
        !==================
        ! Arguments
        class(BoundaryFunctionConstraintsUDT)       :: constraints
        real(R8), intent(in)            :: lambda(:), values(:)
        type(MySparseUDT)               :: hessG, jacG
        type(GridUDT)                   :: grid
        type(MagneticFieldUDT)          :: magneticField
        type(EnvironmentUDT)            :: environment
        class(DesignVariablesGDUDT)     :: designvariables

        ! Auxiliary
        real(R8), allocatable           :: val(:)

        ! Initialize
        !===========
        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            plf     => constraints%plf,         &
            tv      => constraints%vert,        &
            ntv     => constraints%nvert,       &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Compute jacG
        !=============
        ! Simply call differentiation of plf w.r.t. polygonset coordinates
        allocate(val(size(tv)))
        call plf%Evaluate(x(tv), y(tv), 0, 0, val, 'polygonsetcoordinates', values, jacG)

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Destructor

    !------------------------------------------------------------------!
    !                              X-POINTS                            !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeXPointConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the x-point constraints. Here, the x-point is
        ! constrained to its initial location (so it is based on the
        ! initial grid, not yet on the magnetic field flux gradient).

        ! Declare variables
        !==================
        ! Arguments
        class(XPointConstraintsUDT)                 :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                   :: designvariables

        ! Auxiliary
        integer(I8)                                 :: nxpind
        integer(I8), allocatable                    :: xpind(:), order(:)

        logical, allocatable                        :: mask(:)

        ! Data

        ! Initialize
        !===========
        ! Check allocation
        if (allocated(constraints%xpind)) then
            deallocate(constraints%xpind)
        end if

        ! Get x-points
        !=============
        ! Get
        call DetermineXPoints(xpind, nxpind, order, grid)

        ! Allocate
        allocate(mask(nxpind))

        ! Check which x-points are not yet constrained
        mask(:) = .false.
        where (monitor%eqvcc(xpind) <= monitor%maxeqvcc(xpind)-2) mask = .true.

        ! Set method
        constraints%meth = options%xpoptions%meth

        ! Allocate and assign
        constraints%nxpind = count(mask)
        allocate(constraints%xpind(count(mask)))
        constraints%xpind = pack(xpind, mask)

        ! Get current§ coordinates, assign
        constraints%locx = grid%vert%x(constraints%xpind)
        constraints%locy = grid%vert%y(constraints%xpind)

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = 2*constraints%nxpind

        ! Monitor
        monitor%eqvcc(constraints%xpind) = &
            monitor%eqvcc(constraints%xpind) + 2

        ! Debugging info
        !===============
        ! Write datafile
        if (options%writedata == 1) then
            call WriteXPointConstraintVertices(grid, constraints%xpind)
        end if

    end subroutine

    ! Evaluation
    subroutine EvaluateXPointConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! Evaluate the X-point constraints (if location based) as:
        !
        !       G(2*i-1) = (x_i - x_i0)**2 + (x_i - x_i0) = 0
        !       G(2*i)   = (y_i - y_i0)**2 + (y_i - y_i0) = 0
        !
        ! where x_i and y_i are the i-th x-point coordinates, x_i0, y_i0
        ! are the desired (and constant) x-point coordinates.

        ! The gradient is computed as follows:
        !
        !       J(2*i-1, xind(i))      = 2*(x_i - x_i0) + 1
        !       J(2*i, xind(i) + nv)   = 2*(y_i - y_i0) + 1
        !
        ! where nv is the number of vertices in the grid, and xind the
        ! vertex ID vector of all x-points considered.

        ! The Hessian of each ith constraint is:
        !
        !       H(2*i-1, xind(i), xind(i)) = 2
        !       H(2*i, xind(i) + nv, xind(i) + nv) = 2
        !       Otherwise zero
        !
        ! Therefore, the multiplication Hjk,i lambda_i is equal to:
        !
        !       Hjk,i lambda_i = 2 * lambda_i

        ! For gradient based constraints, we simply set dpsidx, dpsidy
        ! to zero. Note that we require 3rd order derivatives then for
        ! hessian computation.

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is
        ! returned, but the hessian-vector multiplication with the
        ! vector lambda, which should be of suitable size.

        ! Note 2: the row and column indices for the constraints are
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine.

        ! Note 3: at first, we compute the linearization of the
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e.
        ! transpose) to obtain the gradient.

        ! Declare variables
        !==================
        ! Arguments
        class(XPointConstraintsUDT)         :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: i, ic, ivg, ivh

        ! Auxiliary variables
        real(R8), allocatable               :: xpx(:), xpy(:), &
            dpsidx(:), dpsidy(:)

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if
        if (present(dGdvarin)) then
            dGdvar = dGdvarin
        end if
        if (present(dgradGdvarin)) then
            dgradGdvar = dgradGdvarin
        else
            dgradGdvarin = SpZeros(0, 0)
        end if

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            tv      => constraints%xpind,       &
            ntv     => constraints%nxpind,      &
            locx    => constraints%locx,        &
            locy    => constraints%locy,        &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateBoundaryFunctionConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Precompute some values if necessary
        select case (trim(constraints%meth))

        case ('loc')

            ! Nothing to compute

        case ('grad')

            ! Extract x-point locations
            xpx = grid%vert%x(constraints%xpind)
            xpy = grid%vert%y(constraints%xpind)

            ! Compute local magnetic field
            allocate(dpsidx(ntv), dpsidy(ntv))
            call magneticField%interp%Evaluate(xpx, xpy, 1, 0, dpsidx)
            call magneticField%interp%Evaluate(xpx, xpy, 0, 1, dpsidy)

        case default

            call gdErrorHandler('Unknown method to impose X-point constraints, choose "loc" or "grad"')

        end select

        ! Evaluate
        select case (trim(constraints%meth))

        case ('loc')

            ! Location based
            do i = 1, ntv
                G(2*i-1) = ( x(tv(i)) - locx(i) )**2 + ( x(tv(i)) - locx(i) )
                G(2*i)   = ( y(tv(i)) - locy(i) )**2 + ( y(tv(i)) - locy(i) )
            end do

        case ('grad')

            ! Gradient based
            G(ic+1:ic+ntv) = dpsidx
            ic = ic + ntv
            G(ic+1:ic+ntv) = dpsidy

        end select

        ! Derivatives
        !====================
        ! Initialize
        gradG%nrow = designvariables%nphi
        gradG%ncol = constraints%ncon
        hessG%nrow = designvariables%nphi
        hessG%ncol = designvariables%nphi

        ! Check which derivatives to compute
        select case (trim(designvariables%type))

        case ('coordinates', 'coordinates_desiredflux')

            ! Only coordinate contributions, no flux contribution
            call constraints%EvaluateCoordinatesDerivative(gradG, hessG, &
                dogradient, dohessian, lambda, grid, magneticField, xpx, xpy)

        case default

            call gdErrorHandler('Unknown design variable type in X-point constraint evaluation')

        end select

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    ! Derivative evaluation
    subroutine EvaluateCoordinatesDerivativeXPointConstraints(&
        constraints, gradG, hessG, dogradient, dohessian, lambda, &
        grid, magneticField, xpx, xpy)

        ! Description
        !============
        ! Evaluate derivatives w.r.t. coordinates (design variable
        ! indices are 'local', meaning we assume the coordinates are
        ! the only design variables)

        ! Declare variables
        !==================
        ! Arguments
        class(XPointConstraintsUDT)         :: constraints
        type(MySparseUDT)                   :: gradG, hessG, jacG
        logical, intent(in)                 :: dogradient, dohessian
        real(R8), allocatable, intent(in)   :: lambda(:)
        type(MagneticFieldUDT)              :: magneticField
        type(GridUDT), intent(in)           :: grid
        real(R8), intent(in)                :: xpx(:), xpy(:)

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               ::  valxx(:), valxy(:), &
            valyy(:), d2psidx2(:), d2psidxdy(:), &
            d2psidy2(:), d3psidx3(:), d3psidx2dy(:), d3psidxdy2(:), &
            d3psidy3(:)

        ! Initialize
        !===========
        associate(&
            nc      => constraints%ncon,        &
            tv      => constraints%xpind,       &
            ntv     => constraints%nxpind,      &
            locx    => constraints%locx,        &
            locy    => constraints%locy,         &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Counters
        ivg = 0
        ivh = 0
        ic = 0

        ! Check case
        !===========
        select case (trim(constraints%meth))

        case ('loc')

            ! Gradient
            if (dogradient) then

                ! Allocate
                jacG%nval = 2*ntv
                call jacG%Allocate()
                allocate(conindex(ntv))
                allocate(valindex(ntv))

                ! x-contribution
                !---------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+2*ntv-1, 2)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv
                jacG%val(valindex) = 2*( x(tv) - locx ) + 1

                ! y-contribution
                !---------------
                ! Build indices
                ivg = ivg + ntv
                conindex = [(k, k = ic+2, ic+2*ntv, 2)]
                valindex = valindex + ntv

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv + grid%vert%ntot
                jacG%val(valindex) = 2*( y(tv) - locy ) + 1

                ! Build gradient
                gradG%nrow = jacG%ncol
                gradG%ncol = jacG%nrow
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            end if

            ! Hessian
            if (dohessian) then

                ! Allocate
                hessG%nval = 4*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(ntv))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))

                ! Compute contributions
                valxx(:) = 2
                valyy(:) = 2
                valxy(:) = 0

                ! xx-contribution
                !----------------
                k = 1
                ! Build indices
                conindex = [(k, k = ic+1, ic+2*ntv-1, 2)]
                valindex = [(k, k = ivh+1, ivh+ntv)]

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv
                hessG%val(valindex) = valxx*lambda(conindex) ! element-wise mult.

                ! xy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = 0 ! element-wise mult.

                ! yx-contribution
                !----------------
                ! symmetric with xy
                ! Build indices
                ivh = ivh + ntv
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv
                hessG%val(valindex) = 0 ! element-wise mult.

                ! yy-contribution
                !----------------
                ! Build indices
                ivh = ivh + ntv
                conindex = [(k, k = ic+2, ic+2*ntv, 2)]
                valindex = valindex + ntv

                ! Add values
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = valyy*lambda(conindex) ! element-wise mult.

            end if

        case ('grad')

            ! Gradient
            if (dogradient) then

                ! Allocate
                jacG%nval = 4*ntv
                call jacG%Allocate()
                allocate(conindex(ntv))
                allocate(valindex(ntv))
                allocate(d2psidx2(ntv), d2psidxdy(ntv), d2psidy2(ntv))

                ! Evaluate
                call magneticField%interp%Evaluate(xpx, xpy, 2, 0, d2psidx2)
                call magneticField%interp%Evaluate(xpx, xpy, 1, 1, d2psidxdy)
                call magneticField%interp%Evaluate(xpx, xpy, 0, 2, d2psidy2)

                ! x-constraints
                !--------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+ntv)]
                valindex = [(k, k = ivg+1, ivg+ntv)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv
                jacG%val(valindex) = d2psidx2
                ivg = ivg + ntv

                valindex = valindex + ntv
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv + grid%vert%ntot
                jacG%val(valindex) = d2psidxdy
                ivg = ivg + ntv

                ! Update counter
                ic = ic + ntv

                ! y-constraints
                !--------------
                ! Build indices
                conindex = [(k, k = ic+1, ic+ntv)]

                ! Add values
                valindex = valindex + ntv
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv
                jacG%val(valindex) = d2psidxdy
                ivg = ivg + ntv

                valindex = valindex + ntv
                jacG%row(valindex) = conindex
                jacG%col(valindex) = tv + grid%vert%ntot
                jacG%val(valindex) = d2psidy2
                ivg = ivg + ntv

                ! Build gradient
                !===============
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            end if

            ! Hessian
            ic = 0
            if (dohessian) then

                ! Allocate
                hessG%nval = 8*ntv
                if (.not. allocated(valindex)) then
                    allocate(valindex(ntv))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(ntv))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(ntv))
                allocate(valxy(ntv))
                allocate(valyy(ntv))
                allocate(d3psidx3(ntv), d3psidx2dy(ntv), d3psidxdy2(ntv), &
                    d3psidy3(ntv))

                ! Compute
                call magneticField%interp%Evaluate(xpx, xpy, 3, 0, d3psidx3)
                call magneticField%interp%Evaluate(xpx, xpy, 2, 1, d3psidx2dy)
                call magneticField%interp%Evaluate(xpx, xpy, 1, 2, d3psidxdy2)
                call magneticField%interp%Evaluate(xpx, xpy, 0, 3, d3psidy3)

                ! x-constraint
                !-------------
                ! Set index
                conindex = [(k, k = ic+1, ic+ntv)]

                ! xx
                valindex = [(k, k = ivh+1, ivh+ntv)]
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv
                hessG%val(valindex) = d3psidx3*lambda(conindex)
                ivh = ivh + ntv

                ! xy
                valindex = valindex + ntv
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = d3psidx2dy*lambda(conindex)
                ivh = ivh + ntv

                ! yx
                valindex = valindex + ntv
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv
                hessG%val(valindex) = d3psidx2dy*lambda(conindex)
                ivh = ivh + ntv

                ! yy
                valindex = valindex + ntv
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = d3psidxdy2*lambda(conindex)
                ivh = ivh + ntv

                ! Update counter
                ic = ic + ntv

                ! y-constraint
                !-------------
                ! Set index
                conindex = [(k, k = ic+1, ic+ntv)]

                ! xx
                valindex = [(k, k = ivh+1, ivh+ntv)]
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv
                hessG%val(valindex) = d3psidx2dy*lambda(conindex)
                ivh = ivh + ntv

                ! xy
                valindex = valindex + ntv
                hessG%row(valindex) = tv
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = d3psidxdy2*lambda(conindex)
                ivh = ivh + ntv

                ! yx
                valindex = valindex + ntv
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv
                hessG%val(valindex) = d3psidxdy2*lambda(conindex)
                ivh = ivh + ntv

                ! yy
                valindex = valindex + ntv
                hessG%row(valindex) = tv + grid%vert%ntot
                hessG%col(valindex) = tv + grid%vert%ntot
                hessG%val(valindex) = d3psidy3*lambda(conindex)
                ivh = ivh + ntv

            end if

        case default

            call gdErrorHandler('Unknown method for imposing X-point constraints')

        end select

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                           EDGE LENGTHS                           !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeEdgeLengthsConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the (desired) edge length constraints. There are
        ! two main parameters here:
        !
        ! - the edge (or face in 2D) indices to be constrained
        ! - the desired length (in [m]) of these edges.
        !
        ! For the former, some preset options are implemented through
        ! routines in the gdmod_utility_optimization module, which can
        ! be called here. These currently include:
        !
        ! - edges that have one boundary vertex that lies on the vessel
        !   (these are basically the edges that determine the width of
        !   the boundary cell near the target plates or other vessel
        !   segments)
        ! - edges that have an x-point as vertex.
        !
        ! For the latter (i.e. the length), it is in principle possible
        ! to apply different lengths for all edges, though the options
        ! here allow only to specify a uniform edge length automatically
        ! which is sufficient for most purposes. Otherwise, this should
        ! be implemented/given manually.

        ! Notes
        !======
        ! Note 1: in determining to which vertex to assing this
        ! constraint (for the monitor), priority is given to the
        ! boundary vertices, if it is possible to assign it there.

        ! Declare variables
        !==================
        ! Arguments
        class(EdgeLengthsConstraintsUDT)            :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                   :: designvariables

        ! Auxiliary
        real(R8)                    :: edgedistvessel, edgedistxpoint

        real(R8), allocatable       :: dvesseledges(:), dxpointedges(:), &
            dgoatedges(:), dx(:), dy(:)

        integer(I8)                 :: nvesseledges, nxpointedges, cc, &
            ev(1:2), ngoatedges

        integer(I8), allocatable    :: vesseledges(:, :), &
            xpointedges(:, :), tempvesseledges(:, :), &
            tempxpointedges(:, :), goatedges(:, :), goatedgeind(:)

        logical                     :: dovesseledges, doxpointedges, &
            doTP, doWG

        logical, allocatable        :: dovesseledgecon(:), &
            doxpointedgecon(:), dogoatedgecon(:)

        ! Loop
        integer(I8)                 :: i, j, k

        ! Data

        ! Initialize
        !===========
        ! Check allocation
        if (allocated(constraints%edgevert)) then
            deallocate(constraints%edgevert)
        end if
        if (allocated(constraints%d)) then
            deallocate(constraints%d)
        end if

        ! Associate
        associate(&
            vert        => grid%vert,           &
            face        => grid%face,           &
            vcc         => monitor%eqvcc,       &
            maxvcc      => monitor%maxeqvcc,    &
            opt         => options%eloptions)

        ! Check if goat data is available
        if (opt%usegoatdata .and. .not. grid%data%hasGoatGGData) then
            print *, 'InitializeEdgeLengthsConstraints: goat data not available, ' // &
                'using user-specified data'
        end if

        ! Do vessel edge lengths?
        dovesseledges   = (opt%dovesseledges == 1)
        edgedistvessel  = opt%edgedistvessel ! desired edge length in [m]
        doTP            = (opt%doTP == 1) ! do target plates?
        doWG            = (opt%doWG == 1) ! do wide grid boundaries?

        ! Do x-point edge lengths?
        doxpointedges   = (opt%doxpointedges == 1)
        edgedistxpoint  = opt%edgedistxpoint ! desired edge length in [m]

        ! Determine edges
        !================
        ! Edges near vessel
        !------------------
        nvesseledges = 0
        if (dovesseledges) then
            call DetermineFluxAlignedVesselEdges(nvesseledges, &
                tempvesseledges, grid, doTP, doWG)

            ! Check for these edges whether they can be constrained
            allocate(dovesseledgecon(nvesseledges))
            dovesseledgecon(:) = .true.
            where ( (vcc(tempvesseledges(:, 1)) >= maxvcc(tempvesseledges(:, 1))) .and. &
                (vcc(tempvesseledges(:, 2)) >= maxvcc(tempvesseledges(:, 2)))) &
                dovesseledgecon = .false.

            ! Skip edges that have two boundary vertices
            where (vert%BV(tempvesseledges(:, 1)) .and. vert%BV(tempvesseledges(:, 2))) dovesseledgecon = .false.

            ! Recompute edges
            nvesseledges = count(dovesseledgecon)
            allocate(vesseledges(nvesseledges, 2))
            vesseledges(:, 1) = pack(tempvesseledges(:, 1), &
                dovesseledgecon)
            vesseledges(:, 2) = pack(tempvesseledges(:, 2), &
                dovesseledgecon)

            ! Set lengths
            allocate(dvesseledges(nvesseledges))
            dvesseledges(:) = edgedistvessel
        end if

        ! Edges near x-point(s)
        !----------------------
        nxpointedges = 0
        if (doxpointedges) then
            call DetermineFluxAlignedXPointEdges(nxpointedges, &
                tempxpointedges, grid)

            ! Check for these edges whether they can be constrained
            allocate(doxpointedgecon(nxpointedges))
            doxpointedgecon(:) = .true.
            where ( (vcc(tempxpointedges(:, 1)) >= maxvcc(tempxpointedges(:, 1))) .and. &
                (vcc(tempxpointedges(:, 2)) >= maxvcc(tempxpointedges(:, 2)))) &
                doxpointedgecon = .false.

            ! Recompute edges
            nxpointedges = count(doxpointedgecon)
            allocate(xpointedges(nxpointedges, 2))
            xpointedges(:, 1) = pack(tempxpointedges(:, 1), &
                doxpointedgecon)
            xpointedges(:, 2) = pack(tempxpointedges(:, 2), &
                doxpointedgecon)

            ! Set lengths
            allocate(dxpointedges(nxpointedges))
            dxpointedges(:) = edgedistxpoint
        end if

        ! Goat edges
        !-----------
        ngoatedges = 0
        if (grid%data%hasGoatGGData .and. opt%usegoatdata) then
            ! Check for these edges whether they can be constrained
            allocate(dogoatedgecon(face%ntot))
            dogoatedgecon = .false.
            where (grid%data%goatggdata%BLind > 0) dogoatedgecon = .true.

            ! Check if we apply current length or user-specified lengths.
            ! For the latter case, we need to check which faces are
            ! not considered
            if (.not. opt%useoriginallength) then
                where (grid%data%goatggdata%BLind > size(opt%edgedistgoat)) dogoatedgecon = .false.
            end if

            ! Don't consider edges that are already fully constrained
            where ( (vcc(face%vert(:, 1)) >= maxvcc(face%vert(:, 1))) .and. &
                (vcc(face%vert(:, 2)) >= maxvcc(face%vert(:, 2)))) &
                dogoatedgecon = .false.

            ! Recompute edges
            ngoatedges = count(dogoatedgecon)
            allocate(goatedges(ngoatedges, 2), goatedgeind(ngoatedges))
            goatedges(:, 1) = pack(face%vert(:, 1), &
                dogoatedgecon)
            goatedges(:, 2) = pack(face%vert(:, 2), &
                dogoatedgecon)
            goatedgeind = pack([(k, k = 1, face%ntot)], dogoatedgecon)

            ! Set lengths
            allocate(dgoatedges(ngoatedges))
            if (opt%useoriginallength) then
                dx = vert%x(goatedges(:, 2)) - vert%x(goatedges(:, 1))
                dy = vert%y(goatedges(:, 2)) - vert%y(goatedges(:, 1))
                dgoatedges = sqrt(dx**2 + dy**2)
            else
                do i = 1, size(opt%edgedistgoat)
                    where (grid%data%goatggdata%BLind(goatedgeind) == i)
                        dgoatedges = opt%edgedistgoat(i)
                    end where
                end do
            end if
        end if

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = nvesseledges + nxpointedges + ngoatedges
        constraints%nedges = nvesseledges + nxpointedges + ngoatedges
        allocate(constraints%edgevert(constraints%nedges, 2))
        allocate(constraints%d(constraints%nedges))

        cc = 0 ! constraint counter
        if (dovesseledges) then
            do j = 1, 2
                constraints%edgevert(cc+1:cc+nvesseledges, j) = &
                    vesseledges(:, j)
            end do
            constraints%d(cc+1:cc+nvesseledges) = &
                dvesseledges
            cc = cc + nvesseledges
        end if
        if (doxpointedges) then
            do j = 1, 2
                constraints%edgevert(cc+1:cc+nxpointedges, j) = &
                    xpointedges(:, j)
            end do
            constraints%d(cc+1:cc+nxpointedges) = &
                dxpointedges
            cc = cc + nxpointedges
        end if
        if (grid%data%hasGoatGGData .and. opt%usegoatdata) then
            do j = 1, 2
                constraints%edgevert(cc+1:cc+ngoatedges, j) = &
                    goatedges(:, j)
            end do
            constraints%d(cc+1:cc+ngoatedges) = &
                dgoatedges
            cc = cc + ngoatedges
        end if

        ! If we use the original length, overwrite
        if (opt%useoriginallength) then
            dx = vert%x(constraints%edgevert(:, 2)) - vert%x(constraints%edgevert(:, 1))
            dy = vert%y(constraints%edgevert(:, 2)) - vert%y(constraints%edgevert(:, 1))
            constraints%d = sqrt(dx**2 + dy**2)
        end if

        ! Monitor
        do i = 1, constraints%nedges
            ! Unpack for ease
            ev = constraints%edgevert(i, 1:2)

            ! Update counter
            if (vert%BV(ev(1)) .and. (vcc(ev(1)) < maxvcc(ev(1)))) then
                ! Assign to boundary vertex
                vcc(ev(1)) = vcc(ev(1)) + 1
            elseif (vert%BV(ev(2)) .and. (vcc(ev(2)) < maxvcc(ev(2)))) then
                ! Assign to boundary vertex
                vcc(ev(2)) = vcc(ev(2)) + 1
            elseif (vcc(ev(1)) <= maxvcc(ev(1))) then
                ! Assign to first vertex, no boundary vertex
                vcc(ev(1)) = vcc(ev(1)) + 1
            elseif (vcc(ev(2)) <= maxvcc(ev(2))) then
                ! Assign to first vertex, no boundary vertex
                vcc(ev(2)) = vcc(ev(2)) + 1
            else
                ! Something wrong - indicates that this edge shouldn't
                ! have been added, though this should've been catched
                ! before
                call gdErrorHandler('InitializeEdgelengthsConstraints:' &
                    // 'edge with vertex IDs as shown above should not' &
                    // ' have been included in the constraints')
            end if
        end do

        ! Debugging info
        !===============
        ! Write datafile
        if (options%writedata == 1) then
            call WriteEdgelengthsConstraintVertexPairs(grid, &
            constraints%edgevert)

        end if

        ! Housekeeping
        !=============
        ! Deallocate
        if (dovesseledges) then
            deallocate(vesseledges, tempvesseledges, dovesseledgecon, &
            dvesseledges)
        end if
        if (doxpointedges) then
            deallocate(xpointedges, tempxpointedges, doxpointedgecon, &
            dxpointedges)
        end if

        ! Deassociate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateEdgeLengthsConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! The edge length constraint for the i-th edge (and therefore
        ! in this case the i-th constraint) can be mathematically
        ! formulated as:
        !
        !       G(i) = 0.5 * (l_i ** 2  - d0 ** 2)
        !
        ! This form ensures that at least one part of the gradient
        ! w.r.t. the coordinate values is non-zero (for non-zero
        ! desired distance, which should always be the case!), while
        ! keeping the constraint formulation as a simple quadratic
        ! function in terms of the vertex coordinates. Here, l_i is
        ! defined as the L2 norm of the edge vector, where the latter is
        ! formed by the vector between its vertices.
        !
        ! The Jacobian w.r.t. the vertex coordinates is then found as
        !
        !       J(i, j) = xv(:, 1) - xv(:, 2), if j = edgevert(i, 1)
        !       J(i, j) = -(xv(:, 1) - xv(:, 2)), if j = edgevert(i, 2)
        !       J(i, j+nv) = yv(:, 1) - yv(:, 2), if j = edgevert(i, 1)
        !       J(i, j+nv) = -(yv(:, 1) - yv(:, 2)), if j = edgevert(i, 2)
        !
        ! where i is the constraint index, j is the vertex index of the
        ! i-th edge (first or second, as specified above), and nv is the
        ! number of grid vertices.
        !
        ! The Hessian w.r.t. the grid vertex coordinates is found as
        !
        !       H(i, j, k) = 1 (j = k)
        !       H(i, j, k) = -1 (j != k)
        !       H(i, j+nv, k+nv) = 1 (j = k)
        !       H(i, j+nv, k+nv) = -1 (j != k)
        !       otherwise zero

        ! Notes
        !======
        ! Note 1: not the true hessian of the constraint vector is
        ! returned, but the hessian-vector multiplication with the
        ! vector lambda, which should be of suitable size.

        ! Note 2: the row and column indices for the constraints are
        ! local, meaning that in no way other constraints are accounted
        ! for in positioning the elements in the matrix. This should be
        ! done in an overarching routine.

        ! Note 3: at first, we compute the linearization of the
        ! constraints, meaning that we actually compute the Jacobian.
        ! Afterwards, we switch the row and column indices (i.e.
        ! transpose) to obtain the gradient.

        ! Declare variables
        !==================
        ! Arguments
        class(EdgeLengthsConstraintsUDT)    :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:)

        ! Auxiliary variables
        real(R8), allocatable               :: valxx(:), valxy(:), &
            valyy(:), xv1(:), xv2(:), yv1(:), yv2(:), dist(:), Gv(:)

        ! Initialize
        !===========
        ! Initialize
        allocate(Gv(grid%vert%ntot))
        Gv = 0.0_R8

        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if
        if (present(dGdvarin)) then
            dGdvar = dGdvarin
        end if
        if (present(dgradGdvarin)) then
            dgradGdvar = dgradGdvarin
        end if

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            nc      => constraints%ncon,        &
            ev      => constraints%edgevert,    &
            d0      => constraints%d,           &
            nv      => grid%vert%ntot,          &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateEdgelengthsConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Compute lengths
        allocate(xv1(nc), xv2(nc), yv1(nc), yv2(nc), dist(nc))
        xv1 = x(ev(:, 1))
        yv1 = y(ev(:, 1))
        xv2 = x(ev(:, 2))
        yv2 = y(ev(:, 2))

        ! Evaluate
        dist = (xv1 - xv2)**2 + (yv1 - yv2)**2
        G = 0.5*(dist - d0**2)
        !do ic = 1, nc
        !    Gv(ev(ic, 1)) = Gv(ev(ic, 1)) + abs(G(ic))
        !    Gv(ev(ic, 2)) = Gv(ev(ic, 2)) + abs(G(ic))
        !end do

        ! Write
        !call Write3DCoordinateData(grid%vert%x, grid%vert%y, Gv, 'con_el_val_vertices')

        ! Constraint gradient
        !====================
        if (dogradient) then
            ! Initialize
            jacG%nrow = nc
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux') ! no flux contributions anyway

                ! Order in jacobian: first x, then y.

                ! Allocate
                jacG%nval = 4*nc
                call jacG%Allocate()
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! x-contribution
                !---------------
                ! Build indices for xv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 1)
                jacG%val(valindex) = xv1 - xv2

                ! Update counters
                ivg = ivg + nc

                ! Build indices for xv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 2)
                jacG%val(valindex) = -(xv1 - xv2)

                ! Update counters
                ivg = ivg + nc

                ! y-contribution
                !---------------
                ! Build indices for yv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 1) + grid%vert%ntot
                jacG%val(valindex) = yv1 - yv2

                ! Update counters
                ivg = ivg + nc

                ! Build indices for yv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 2) + grid%vert%ntot
                jacG%val(valindex) = -(yv1 - yv2)

                ! Update counters
                ivg = ivg + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol
                gradG%ncol = jacG%nrow
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi
            hessG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux')

                ! Allocate
                hessG%nval = 8*nc
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(nc))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(nc))
                allocate(valxy(nc))
                allocate(valyy(nc))

                ! Compute contributions
                valxx(:) = 1
                valyy(:) = 1
                valxy(:) = 0

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! xx-contribution
                !----------------
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 1)
                hessG%col(valindex) = ev(:, 1)
                hessG%val(valindex) = valxx*lambda(conindex) ! x1x1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 1)
                hessG%col(valindex) = ev(:, 2)
                hessG%val(valindex) = -valxx*lambda(conindex) ! x1x2

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 2)
                hessG%col(valindex) = ev(:, 1)
                hessG%val(valindex) = -valxx*lambda(conindex) ! x2x1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 2)
                hessG%col(valindex) = ev(:, 2)
                hessG%val(valindex) = valxx*lambda(conindex) ! x2x2

                ! Update counters
                ivh = ivh + nc

                ! xy-contribution
                !----------------
                ! no contributions

                ! yx-contribution
                !----------------
                ! no contributions

                ! yy-contribution
                !----------------
                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 1) + nv
                hessG%col(valindex) = ev(:, 1) + nv
                hessG%val(valindex) = valxx*lambda(conindex) ! y1y1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 1) + nv
                hessG%col(valindex) = ev(:, 2) + nv
                hessG%val(valindex) = -valxx*lambda(conindex) ! y1y2

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 2) + nv
                hessG%col(valindex) = ev(:, 1) + nv
                hessG%val(valindex) = -valxx*lambda(conindex) ! y2y1

                ! Update counters
                ivh = ivh + nc

                ! Build indices
                valindex = [(k, k = ivh+1, ivh+nc)]

                ! Add values
                hessG%row(valindex) = ev(:, 2) + nv
                hessG%col(valindex) = ev(:, 2) + nv
                hessG%val(valindex) = valxx*lambda(conindex) ! y2y2

                ! Update counters
                ivh = ivh + nc

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate
        if (dogradient) then
            deallocate(valindex, conindex)
        end if

        if (dohessian) then
            if (allocated(valindex)) then
                deallocate(valindex, conindex)
            end if
            deallocate(valxx, valxy, valyy)
        end if

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                           ORTHOGONALITY                          !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeOrthogonalityConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the edge length constraints on the desired edges.
        ! As there are many ways to determine the edges, this is
        ! cast into a function located in the gdmod_utility_optimization
        ! module, i.e. DetermineEdgesOrthogonalityConstraints. Different
        ! preset options are present:
        !
        ! - simple boxes that determine, starting from the original
        !   grid, which edges to be constrained (or not) by checking if
        !   edges are within the box.
        ! - flux value based: edges with a flux value between certain
        !   limits are constrained (useful for e.g. the core)
        ! - initial orthogonality based: numerically compute the
        !   deviation from orthogonality of the original edges and
        !   decide based on that whether to include the edges.
        !
        ! See the routine for more options and details.

        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(OrthogonalityConstraintsUDT)          :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                   :: designvariables

        ! Auxiliary
        integer(I8)                 :: nincludebox, nexcludebox, tv, &
            startindex, endindex, tID, nbID(1:2), vpc, tbv, tnbv, tf
        real(R8)                    :: tx, ty, tn, bx, by, bn, dotprod, &
            epsperp, maxx, minx, maxy, miny
        logical                     :: checkperp, isfaceperp, &
            debugplots, delv(1:2)

        integer(I8), allocatable    :: cvertlist(:), temp(:), &
            northcon(:), maxnorthcon(:), tvn(:), vpairs(:, :), &
            tvn1(:), tvn2(:), tvn_old(:), closuretype(:)

        real(R8), allocatable       :: Btx(:), Bty(:), xf(:), yf(:)

        logical, allocatable        :: cvert(:), boxcheck(:), &
            movetoback(:), movetofront(:), ismarked(:), isconstrained(:), &
            isperp(:), cID(:), isvesselvertex(:), isvesselface(:)

        ! Loop
        integer(I8)                 :: i, j, k

        ! Data

        ! Initialize
        !===========
        ! Check allocation
        if (allocated(constraints%edgevert)) then
            deallocate(constraints%edgevert)
        end if

        ! Associate
        associate(&
            interp      => magneticField%interp,    &
            vert        => grid%vert,               &
            vcc         => monitor%eqvcc,           &
            maxvcc      => monitor%maxeqvcc,        &
            x           => grid%vert%x,             &
            y           => grid%vert%y,             &
            opt         => options%orthoptions)

        ! Debug plots?
        debugplots = .false.

        ! Boxes for edges to be included?
        nincludebox     = size(opt%includeboxx, 1) ! number of boxes

        ! Boxes for edges to be excluded
        nexcludebox     = size(opt%excludeboxx, 1) ! number of boxes

        ! Tolerances
        epsperp     = opt%epsperp

        ! Check perpendicularity?
        checkperp = (opt%checkperp == 1)

        ! Get vessel vertices
        call DetermineVesselVertices(isvesselvertex, isvesselface, grid)

        ! Determine edge vertices
        !========================
        ! Compute magnetic field vector components Btx, Bty
        allocate(Btx(vert%ntot), Bty(vert%ntot))
        call interp%Evaluate(x, y, 0, 1, Btx)
        call interp%Evaluate(x, y, 1, 0, Bty)
        Btx = -Btx ! adjust sign

        ! Determine which nodes to consider
        allocate(cvert(vert%ntot))
        allocate(boxcheck(vert%ntot))
        allocate(northcon(vert%ntot), maxnorthcon(vert%ntot))

        cvert(:) = .false.
        northcon(:) = 0
        maxnorthcon(:) = 2
        where (vert%BV) maxnorthcon = 1

        ! Include boundary vertices with zero ID?
        if (opt%includecutcellvert) then
            where (isvesselvertex) cvert = .true.
        end if

        ! Include core vertices?
        if (opt%includecorevert) then
            ! Check which flux surfaces are closed
            closuretype = DetermineFluxSurfaceClosure(grid)
            do i = 1, size(closuretype)
                if (closuretype(i) == 1) then
                    ! Closed flux surface
                    where (grid%vert%fieldlineID == i) cvert = .true.
                end if
            end do
        end if

        ! Include?
        do i = 1, nincludebox ! include points in the box
            maxx = maxval(opt%includeboxx(i, :))
            minx = minval(opt%includeboxx(i, :))
            maxy = maxval(opt%includeboxy(i, :))
            miny = minval(opt%includeboxy(i, :))
            boxcheck = IsInBox(minx, maxx, miny, maxy, vert%x, vert%y)
            cvert = cvert .or. boxcheck
        end do

        ! Exclude?
        do i = 1, nexcludebox ! exclude points outside the box
            maxx = maxval(opt%excludeboxx(i, :))
            minx = minval(opt%excludeboxx(i, :))
            maxy = maxval(opt%excludeboxy(i, :))
            miny = minval(opt%excludeboxy(i, :))
            boxcheck = IsInBox(minx, maxx, miny, maxy, vert%x, vert%y)
            where (boxcheck) cvert = .false.
        end do

        ! Prioritize inner vertices (typically yields better results,
        ! but may not be a general approach)
        allocate(cvertlist(count(cvert)))
        allocate(ismarked(vert%ntot))
        allocate(movetoback(count(cvert)))
        allocate(movetofront(count(cvert)))

        ismarked(:) = .false.
        where (cvert) ismarked = .true.
        movetoback(:) = .false.
        movetofront(:) = .false.
        cvertlist = pack( [(k, k = 1, vert%ntot)], cvert) ! node indices

        do i = 1, count(cvert)
            ! Get the current vertex
            tv = cvertlist(i)

            ! Get the neighbours
            allocate(tvn(vert%neigP(tv, 2)))
            tvn = vert%neig(vert%neigP(tv, 1):vert%neigP(tv, 1)+vert%neigP(tv, 2)-1)

            ! Check if any are not marked
            if (any(.not. ismarked(tvn))) then
                movetoback(i) = .true.
            end if

            ! Check if it is an x-point
            if (.not. movetoback(i)) then
                if (vert%neigP(tv, 2) > 4) then ! crude check for x-point
                    ! Prioritize, move to front
                    movetofront(i) = .true.
                end if
            end if

            ! Deallocate
            deallocate(tvn)
        end do

        ! Rebuild cvertlist
        allocate(temp(size(cvertlist, 1)))
        temp = cvertlist
        startindex = 1
        endindex = count(movetofront)
        cvertlist(startindex:endindex) = pack(temp, movetofront)
        startindex = startindex + endindex
        endindex = endindex + count((.not. movetoback) .and. (.not. movetofront))
        cvertlist(startindex:endindex) = pack(temp, (.not. movetoback) .and. (.not. movetofront))
        startindex = endindex + 1
        endindex = endindex + count(movetoback)
        cvertlist(startindex:endindex) = pack(temp, movetoback)

        deallocate(temp)
        deallocate(movetoback, movetofront)

        ! Initialize
        allocate(isconstrained(grid%face%ntot))
        isconstrained(:) = .false.

        ! Loop
        vpc = 0 ! face counter
        allocate(vpairs(grid%face%ntot, 2))! allocate too big, trim later

        do i = 1, size(cvertlist, 1)
            ! Get the current vertex
            tv = cvertlist(i)

            ! Get the neighbours of this vertex
            allocate(tvn(vert%neigP(tv, 2)))
            tvn = vert%neig(vert%neigP(tv, 1):vert%neigP(tv, 1)+vert%neigP(tv, 2)-1)

            ! Get the fieldline ID
            tID = vert%fieldlineID(tv)

            ! Check if zero
            if (tID .eq. 0) then ! this is a vertex without fieldline
                ! Check which faces are aligned
                allocate(isperp(size(tvn, 1)))

                ! Loop over all vertex neighbours
                do j = 1, size(tvn, 1)
                    ! Get normalized face vector
                    tx = x(tvn(j)) - x(tv)
                    ty = y(tvn(j)) - y(tv)
                    tn = (tx**2 + tx**2)**0.5
                    tx = tx/tn
                    ty = ty/tn

                    ! Get normalized magnetic field vector
                    bx = (Btx(tv) + Btx(tvn(j)))*0.5
                    by = (Bty(tv) + Bty(tvn(j)))*0.5
                    bn = (bx**2 + by**2)**0.5
                    bx = bx/bn ! normalize
                    by = by/bn

                    ! Compute dot product
                    dotprod = tx*bx + ty*by

                    ! Check
                    if (abs(dotprod) < epsperp) then
                        isperp(j) = .true.
                    endif

                end do

                ! Retain perpendicular faces only if 2 found
                if (count(isperp) == 2) then
                    allocate(temp(size(tvn, 1)))
                    temp = tvn
                    deallocate(tvn)
                    allocate(tvn(count(isperp)))
                    tvn = pack(temp, isperp)

                    deallocate(temp)
                else
                    deallocate(tvn)
                    allocate(tvn(0))
                end if

                ! Deallocate
                deallocate(isperp)
            else
                ! Check which vertices have the same ID
                allocate(cID(size(tvn, 1)))
                cID = tID .eq. vert%fieldlineID(tvn)

                ! Extract vertices that have NOT the same ID
                allocate(temp(size(tvn, 1)))
                temp = tvn
                deallocate(tvn)
                allocate(tvn(count(.not. cID)))
                tvn = pack(temp, .not. cID)
                deallocate(temp)

                ! If two nodes remain, and they have different field
                ! line IDs, then constrain. Otherwise don't (probably
                ! stacked triangles)
                if (size(tvn, 1) == 2) then
                    nbID = vert%fieldlineID(tvn)
                    if (nbID(1) == nbID(2)) then
                        deallocate(tvn)
                        allocate(tvn(0))
                    else
                        ! Check if the neighbouring vertices only have
                        ! one vertex with ID equal to the ID of the first
                        ! vertex
                        tvn1 = GetVertNeig(vert, tvn(1))
                        tvn2 = GetVertNeig(vert, tvn(2))
                        delv = .false.
                        if (count(vert%fieldlineID(tvn1) == tID) /= 1) then
                            ! Don't include
                            delv(1) = .true.
                        end if
                        if (isvesselvertex(tvn(1))) then
                            delv(1) = .true.
                        end if
                        if (count(vert%fieldlineID(tvn2) == tID) /= 1) then
                            ! Don't include
                            delv(2) = .true.
                        end if
                        if (isvesselvertex(tvn(2))) then
                            delv(2) = .true.
                        end if
                        tvn_old = tvn
                        deallocate(tvn)
                        allocate(tvn(count(.not. delv)))
                        tvn = pack(tvn_old, .not. delv)
                    end if
                else
                    ! Check if the current vertex is a vessel vertex
                    if (isvesselvertex(tv)) then
                        ! Eliminate vessel vertex neighbours
                        temp = tvn
                        deallocate(tvn)
                        allocate(tvn(count(.not. isvesselvertex(temp))))
                        tvn = pack(temp, .not. isvesselvertex(temp))
                        deallocate(temp)

                        ! Check if one remains, otherwise deallocate
                        if (size(tvn) /= 1) then
                            deallocate(tvn)
                            allocate(tvn(0))
                        end if
                    else
                        deallocate(tvn)
                        allocate(tvn(0))
                    end if
                end if

                ! Deallocate
                deallocate(cID)
            end if

            ! Don't put orthogonality constraints for boundary vertices
            ! with non-zero flux surface - may lead to more unstable
            ! problem and is typically not required/necessary for grid quality
            if (isvesselvertex(tv) .and. (tID /= 0)) then
               deallocate(tvn)
                allocate(tvn(0))
            end if

            ! Constrain each pair
            do j = 1, size(tvn, 1)
                ! Get the face
                call MapVertexPairToFace(tv, tvn(j), grid%face%vert, &
                    grid%face%ntot, tf)

                ! Check initial perpendicularity
                isfaceperp = .true.
                if (checkperp) then
                    ! Get normalized face vector
                    tx = x(tvn(j)) - x(tv)
                    ty = y(tvn(j)) - y(tv)
                    tn = (tx**2 + tx**2)**0.5
                    tx = tx/tn
                    ty = ty/tn

                    ! Get normalized magnetic field vector
                    bx = (Btx(tv) + Btx(tvn(j)))*0.5
                    by = (Bty(tv) + Bty(tvn(j)))*0.5
                    bn = (bx**2 + by**2)**0.5
                    bx = bx/bn ! normalize
                    by = by/bn

                    ! Compute dot product
                    dotprod = tx*bx + ty*by

                    ! Check
                    if (abs(dotprod) > epsperp) then
                        ! Not perpendicular, don't consider this face
                        isfaceperp = .false.
                    end if
                end if

                ! Check if there is a boundary vertex in this edge which
                ! already has been constrained.
                if ( (.not. isconstrained(tf)) & ! should not be constrained
                    .and. ( (vcc(tvn(j)) < maxvcc(tvn(j))) .or. (vcc(tv) < maxvcc(tv)) ) & ! a vertex has less than 2 constraints already imposed
                    .and. ( (.not. vert%BV(tvn(j))) .or. (.not. vert%BV(tv))) & ! at least one is an internal vertex
                    .and. ( isfaceperp ) & ! the face is initially almost orthogonal
                    .and. ( ( northcon(tvn(j)) < maxnorthcon(tvn(j)) ) .and. ( northcon(tv) < maxnorthcon(tv) ) ) &
                    ) then

                    ! Add
                    vpc = vpc + 1
                    vpairs(vpc, :) = [tv, tvn(j)]

                    ! Update counter - first attribute to BV if possible
                    if (vert%BV(tv) .or. vert%BV(tvn(j)) ) then
                        if (vert%BV(tv)) then
                            tbv = tv
                            tnbv = tvn(j)
                        else
                            tbv = tvn(j)
                            tnbv = tv
                        end if

                        ! Check if we can add the constraint there
                        if ( ( northcon(tbv) < maxnorthcon(tbv) ) &
                            .and. vcc(tbv) < maxvcc(tbv) ) then
                            vcc(tbv) = vcc(tbv) + 1
                        else
                            vcc(tnbv) = vcc(tnbv) + 1
                        end if
                    else
                        if ( vcc(tv) < maxvcc(tv) ) then
                            vcc(tv) = vcc(tv) + 1
                        else
                            vcc(tvn(j)) = vcc(tvn(j)) +  1
                        end if
                    end if

                    ! Update counters
                    northcon(tvn(j)) = northcon(tvn(j)) + 1
                    northcon(tv) = northcon(tv) + 1
                    isconstrained(tf) = .true.

                end if

            end do

            ! Deallocate
            deallocate(tvn)
        end do

        ! Update constraint quantities
        !=============================
        ! Constraints
        constraints%ncon = vpc
        constraints%nedges = vpc
        allocate(constraints%edgevert(constraints%nedges, 2))
        constraints%edgevert = vpairs(1:vpc, :)

        ! Construct radial lines
        call constraints%radiallines%Construct(constraints%edgevert, vert%x, vert%y)

        ! Visualize
        if (debugplots) then

            ! Plot the faces that are constrained

            ! Compute face centers
            allocate(xf(vpc), yf(vpc))
            xf = 0.5*(x(constraints%edgevert(:, 1)) + x(constraints%edgevert(:, 2)))
            yf = 0.5*(y(constraints%edgevert(:, 1)) + y(constraints%edgevert(:, 2)))
            call PlotGridWithPoints(grid, xf, yf, '-p')
            deallocate(xf, yf)

        end if

        ! Debugging info
        !===============
        ! Write datafile
        if (options%writedata == 1) then
            call WriteOrthogonalityConstraintVertexPairs(grid, &
                constraints%edgevert)
        end if

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(Btx, Bty, cvert, boxcheck, vpairs)

        ! Deassociate
        end associate

    end subroutine

    ! Evaluation
    subroutine EvaluateOrthogonalityConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! The orthogonality constraints are evaluated per vertex pair,
        ! where one tries to satisfy the following condition:
        !
        !       gx*dx + gy*dy = 0.
        !
        ! This ensures orthogonality of the face w.r.t. the magnetic
        ! field, given by its components gx and gy. dx and dy are the
        ! face tangents (they can be normalized, but this is not
        ! strictly necessary, same for the magnetic field. This does,
        ! however, influence the absolute value of the constraint).


        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(OrthogonalityConstraintsUDT)  :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: i, ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:), &
            row(:), col(:)

        ! Auxiliary variables
        real(R8), allocatable               :: valxx(:), valxy(:), valyx(:), &
            valyy(:), xv(:, :), yv(:, :), xf(:), yf(:), gxf(:), gyf(:), &
            dx(:), dy(:), gxxf(:), gxyf(:), gyxf(:), gyyf(:), gxxxf(:), &
            gxyxf(:), gyxxf(:), gyyxf(:), gxxyf(:), gxyyf(:), gyxyf(:), &
            gyyyf(:), Gv(:)

        ! Initialize
        !===========
        ! Initialize
        allocate(Gv(grid%vert%ntot))
        Gv = 0.0_R8

        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            interp  => magneticField%interp,    &
            nc      => constraints%ncon,        &
            ev      => constraints%edgevert,    &
            vert    => grid%vert,               &
            nv      => grid%vert%ntot,          &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateEdgelengthsConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Precompute
        allocate(xv(nc, 2), yv(nc, 2), xf(nc), yf(nc), gxf(nc), gyf(nc), &
            dx(nc), dy(nc))

        xv(:, 1) = x(ev(:, 1))
        xv(:, 2) = x(ev(:, 2))
        yv(:, 1) = y(ev(:, 1))
        yv(:, 2) = y(ev(:, 2))
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*sum(xv, 2)
        yf = 0.5*sum(yv, 2)
        call interp%Evaluate(xf, yf, 0, 1, gxf)
        call interp%Evaluate(xf, yf, 1, 0, gyf)
        gxf = -gxf ! take correct sign

        ! Evaluate
        G = gxf*dx + gyf*dy
        do i = 1, nc
            Gv(ev(i, 1)) = Gv(ev(i, 1)) + abs(G(i))
            Gv(ev(i, 2)) = Gv(ev(i, 2)) + abs(G(i))
        end do

        ! Write
        !call Write3DCoordinateData(grid%vert%x, grid%vert%y, Gv, 'con_orth_val_vertices')

        ! Constraint gradient
        !====================
        if (dogradient) then
            ! Initialize
            jacG%nrow = nc
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux') ! no flux contributions

                ! Order in jacobian: first x, then y.

                ! Allocate
                jacG%nval = 4*nc
                call jacG%Allocate()
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Precompute
                allocate(gxxf(nc), gxyf(nc), gyxf(nc), gyyf(nc))
                call interp%Evaluate(xf, yf, 1, 1, gxxf)
                call interp%Evaluate(xf, yf, 0, 2, gxyf)
                call interp%Evaluate(xf, yf, 2, 0, gyxf)
                call interp%Evaluate(xf, yf, 1, 1, gyyf)
                gxxf = -gxxf ! correct sign
                gxyf = -gxyf ! correct sign

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! x-contribution
                !---------------
                ! Build indices for xv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 1)
                jacG%val(valindex) = (-gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)

                ! Update counters
                ivg = ivg + nc

                ! Build indices for xv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 2)
                jacG%val(valindex) = (gxf + 0.5*dx*gxxf + 0.5*dy*gyxf)

                ! Update counters
                ivg = ivg + nc

                ! y-contribution
                !---------------
                ! Build indices for yv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 1) + grid%vert%ntot
                jacG%val(valindex) = (-gyf + 0.5*dy*gyyf + 0.5*dx*gxyf)

                ! Update counters
                ivg = ivg + nc

                ! Build indices for yv2
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = ev(:, 2) + grid%vert%ntot
                jacG%val(valindex) = (gyf + 0.5*dy*gyyf + 0.5*dx*gxyf)

                ! Update counters
                ivg = ivg + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol
                gradG%ncol = jacG%nrow
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi
            hessG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux')

                ! Initialize
                !===========
                ! Allocate
                hessG%nval = 16*nc
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(nc))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if
                allocate(valxx(4*nc))
                allocate(valxy(4*nc))
                allocate(valyx(4*nc))
                allocate(valyy(4*nc))
                allocate(row(4*nc), col(4*nc))

                ! Precompute
                allocate(gxxxf(nc), gxyxf(nc), gyxxf(nc), gyyxf(nc), &
                    gxxyf(nc), gxyyf(nc), gyxyf(nc), gyyyf(nc))

                call interp%Evaluate(xf, yf, 2, 1, gxxxf)
                call interp%Evaluate(xf, yf, 1, 2, gxyxf)
                call interp%Evaluate(xf, yf, 3, 0, gyxxf)
                call interp%Evaluate(xf, yf, 2, 1, gyyxf)
                call interp%Evaluate(xf, yf, 1, 2, gxxyf)
                call interp%Evaluate(xf, yf, 0, 3, gxyyf)
                call interp%Evaluate(xf, yf, 2, 1, gyxyf)
                call interp%Evaluate(xf, yf, 1, 2, gyyyf)

                !call EvaluateBicubicSplineInterpolant(xf, yf, gxxxf, &
                !    interp, '2', '1')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gxyxf, &
                !    interp, '1', '2')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gyxxf, &
                !    interp, '2', '1')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gyyxf, &
                !    interp, '3', '0')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gxxyf, &
                !    interp, '1', '2')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gxyyf, &
                !    interp, '0', '3')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gyxyf, &
                !    interp, '1', '2')
                !call EvaluateBicubicSplineInterpolant(xf, yf, gyyyf, &
                !    interp, '2', '1')
                gxxxf = -gxxxf
                gxyxf = -gxyxf
                gxxyf = -gxxyf
                gxyyf = -gxyyf

                ! Compute contributions
                !======================

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! v1, v1
                valindex = [(k, k = ivh+1, ivh+nc)]
                row(valindex) = ev(:, 1)
                col(valindex) = ev(:, 1)
                valxx(valindex) = (0.25*dx*gxxxf - 1.0*gxxf + 0.25*dy*gyxxf)*lambda ! x1x1
                valxy(valindex) = (0.25*dx*gxyxf - 0.5*gyxf - 0.5*gxyf + 0.25*dy*gyyxf)*lambda ! x1y1
                valyx(valindex) = (0.25*dx*gxyxf - 0.5*gyxf - 0.5*gxyf + 0.25*dy*gyyxf)*lambda ! y1x1
                valyy(valindex) = (0.25*dx*gxyyf - 1.0*gyyf + 0.25*dy*gyyyf)*lambda ! y1y1
                ivh = ivh + nc

                ! v1, v2
                valindex = [(k, k = ivh+1, ivh+nc)]
                row(valindex) = ev(:, 1)
                col(valindex) = ev(:, 2)
                valxx(valindex) = (0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x1x2
                valxy(valindex) = (0.5*gyxf - 0.5*gxyf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x1y2
                valyx(valindex) = (0.5*gxyf - 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y1x2
                valyy(valindex) = (0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y1y2
                ivh = ivh + nc

                ! v2, v1
                valindex = [(k, k = ivh+1, ivh+nc)]
                row(valindex) = ev(:, 2)
                col(valindex) = ev(:, 1)
                valxx(valindex) = (0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x2x1
                valxy(valindex) = (0.5*gxyf - 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x2y1
                valyx(valindex) = (0.5*gyxf - 0.5*gxyf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y2x1
                valyy(valindex) = (0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y1y2
                ivh = ivh + nc

                ! v2, v2
                valindex = [(k, k = ivh+1, ivh+nc)]
                row(valindex) = ev(:, 2)
                col(valindex) = ev(:, 2)
                valxx(valindex) = (1.0*gxxf + 0.25*dx*gxxxf + 0.25*dy*gyxxf)*lambda ! x2x2
                valxy(valindex) = (0.5*gxyf + 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! x2y2
                valyx(valindex) = (0.5*gxyf + 0.5*gyxf + 0.25*dx*gxyxf + 0.25*dy*gyyxf)*lambda ! y2x2
                valyy(valindex) = (1.0*gyyf + 0.25*dx*gxyyf + 0.25*dy*gyyyf)*lambda ! y2y2
                ivh = ivh + nc

                ! Build full hessian
                !===================
                hessG%row = [row, row, row+vert%ntot, row+vert%ntot]
                hessG%col = [col, col+vert%ntot, col, col+vert%ntot]
                hessG%val = [valxx, valxy, valyx, valyy]

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Deallocate
        deallocate(gxf, gyf, xv, yv, xf, yf, dx, dy)

        if (dogradient) then
            deallocate(valindex, conindex)
            deallocate(gxxf, gxyf, gyxf, gyyf)
        end if

        if (dohessian) then
            if (allocated(valindex)) then
                deallocate(valindex, conindex)
            end if
            deallocate(valxx, valxy, valyx, valyy)
            deallocate(gxxxf, gxyxf, gyxxf, gyyxf, gxxyf, gxyyf, &
                gyxyf, gyyyf)
        end if

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                         FIXED FLUX VALUES                        !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeFixedFluxvaluesConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the fixed flux value constraints. These constraints
        ! comprise constraints on the core flux and constraints on the
        ! outermost flux surfaces.

        ! Notes
        !======
        ! Note 1: the design variable indices, psiind, are only
        ! allocated here but determined later when finalizing the
        ! optimization problem (see FinalizeInitialization in the
        ! gdmod_optimizationengine module)

        ! Modules
        use mod_definitions, only : coreID, outerboundaryID

        ! Declare variables
        !==================
        ! Arguments
        class(FixedFluxvaluesConstraintsUDT)        :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                 :: designvariables

        ! Auxiliary
        integer(I8)                     :: nfs, ntfsIDs, nfsv, ncfs, &
            ntpind, vID
        integer(I8), allocatable        :: allIDs(:), fsind(:), &
            tvID(:), tfsIDs(:), allvertIDs(:), tvn(:), tv(:), tpind(:), &
            tptype(:), tvnID(:), uniqueID(:)

        logical                         :: dowarning, &
            islonely
        logical, allocatable            :: doesIDoccur(:), &
            hasbeenfound(:), isvesselvertex(:), isvesselface(:), &
            istangencypoint(:)

        real(R8), allocatable           :: psid(:), psi(:)

        ! Loop
        integer(I8)                     :: i, j, k, fscc

        ! Checks
        !=======
        ! Allocation status
        if (allocated(constraints%psiind)) then
            deallocate(constraints%psiind)
        end if

        ! Are the design variable compatible?
        select case (designvariables%type)

        case ('coordinates_desiredflux')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('InitializeFixedFluxvaluesConstraints: design variable type is incompatible')

        end select

        ! Initialize
        !===========
        ! Associate
        associate(&
            maxccv          => monitor%maxeqvcc,        &
            ccv             => monitor%eqvcc,           &
            x               => grid%vert%x,             &
            y               => grid%vert%y,             &
            fieldlineID     => grid%vert%fieldlineID,   &
            docoreflux      => options%ffvoptions%fixcoreflux,  &
            dotpflux        => options%ffvoptions%fixtpflux,    &
            doouterflux     => options%ffvoptions%fixouterflux  &
            )

        ! Allocate
        nfs = maxval(fieldlineID) ! should provide upper bound
        allIDs = [(k, k = 1, nfs)]
        allvertIDs = [(k, k = 1, grid%vert%ntot)]
        fscc = 0 ! flux surface counter
        allocate(fsind(nfs), doesIDoccur(nfs), psid(nfs), &
            hasbeenfound(nfs), psi(grid%vert%ntot))
        hasbeenfound = .false.

        ! Initialize
        fsind = 0
        psid = 0

        ! Get tangency points
        call DetermineTangencyPoints(tpind, ntpind, tptype, grid)
        allocate(istangencypoint(grid%vert%ntot))
        istangencypoint(tpind) = .true.

        ! Core flux
        !==========
        if (docoreflux) then
            ! Determine core flux surface ID(s)
            do i = 1, size(grid%Bnd, 1)
                if (grid%Bnd(i)%ID == coreID) then ! core boundary ID hard coded here...

                    ! Get the fieldline IDs of these vertices
                    tvID = fieldlineID(grid%Bnd(i)%vert)

                    ! Check which ones occur
                    doesIDoccur = .false.
                    doesIDoccur(tvID) = .true.

                    ! Extract IDs
                    ntfsIDs = count(doesIDoccur .and. (.not. hasbeenfound))
                    allocate(tfsIDs(ntfsIDs))
                    tfsIDs = pack(allIDs, doesIDoccur .and. (.not. hasbeenfound))

                    ! Add
                    fsind(fscc+1:fscc+ntfsIDs) = tfsIDs
                    hasbeenfound(tfsIDs) = .true.

                    ! Update counter
                    fscc = fscc + ntfsIDs

                    ! Housekeeping
                    deallocate(tfsIDs)
                end if
            end do

            ! Determine flux value
            select case (options%ffvoptions%fixcorefluxmeth)

            case ('auto')

                ! Precompute flux values
                call magneticField%interp%Evaluate(x, y, 0, 0, psi)

                ! Determine flux value as mean of current flux values
                do i = 1, fscc
                    nfsv = 0
                    do j = 1, grid%vert%ntot
                        if (fieldlineID(j) == fsind(i)) then

                            ! Compute
                            psid(i) = psid(i) + psi(j)
                            nfsv = nfsv + 1

                        end if
                    end do

                    ! Check
                    if (nfsv == 0) then
                        ! No vertices found, throw error
                        print *, 'core flux surface ID: ', fsind(i)
                        call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                            'no vertices found for core flux surface ID')
                    end if

                    ! Average
                    psid(i) = psid(i)/nfsv

                end do

            case ('manual')

                ! Need to specify as much core values are there are
                ! flux surfaces
                if (size(options%ffvoptions%corefluxval, 1) /= fscc) then
                    print *, 'number of core flux surfaces: ', fscc
                    print *, 'field line IDs of core flux surfaces: ', fsind(1:fscc)
                    call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                        'need to specify the amount of core flux surface values mentioned above')
                end if

                ! Add
                psid(1:fscc) = options%ffvoptions%corefluxval

            case default

                ! Throw error
                call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                    'unknown method to determine flux values for core')

            end select

        end if

        ! Store number of core flux surfaces
        ncfs = fscc

        ! Outer flux
        !===========
        if (doouterflux) then
            ! Determine outer flux surface ID(s)
            do i = 1, size(grid%Bnd, 1)
                if (any([grid%Bnd(i)%ID == [outerboundaryID]])) then ! outer boundary ID hard coded here...

                    ! Get the fieldline IDs of these vertices
                    tvID = fieldlineID(grid%Bnd(i)%vert)

                    ! Skip if there is a zero ID - tangency point
                    if (any(tvID == 0_I8)) then
                        cycle
                    end if

                    ! Check which ones occur
                    doesIDoccur = .false.
                    doesIDoccur(tvID) = .true.

                    ! Extract IDs
                    ntfsIDs = count(doesIDoccur .and. (.not. hasbeenfound))
                    allocate(tfsIDs(ntfsIDs))
                    tfsIDs = pack(allIDs, doesIDoccur .and. (.not. hasbeenfound))

                    ! Add
                    fsind(fscc+1:fscc+ntfsIDs) = tfsIDs
                    hasbeenfound(tfsIDs) = .true.

                    ! Update counter
                    fscc = fscc + ntfsIDs

                    ! Housekeeping
                    deallocate(tfsIDs)
                end if
            end do

            ! Determine flux value
            select case (options%ffvoptions%fixouterfluxmeth)

            case ('auto')

                ! Precompute flux values
                call magneticField%interp%Evaluate(x, y, 0, 0, psi)

                ! Precompute vessel vertices
                call DetermineVesselVertices(isvesselvertex, isvesselface, grid)

                ! Determine flux value as mean of current flux values
                do i = ncfs+1, fscc
                    nfsv = 0
                    do j = 1, grid%vert%ntot
                        if ( (fieldlineID(j) == fsind(i)) .and. (isvesselvertex(j))) then

                            ! Compute
                            psid(i) = psid(i) + psi(j)
                            nfsv = nfsv + 1

                        end if
                    end do

                    if (nfsv == 0) then
                        ! No vertices found, throw error
                        print *, 'outer flux surface ID: ', fsind(i)
                        call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                            'no boundary vertices found for outer flux surface ID')
                    end if

                    ! Average
                    psid(i) = psid(i)/nfsv

                end do

            case ('manual')

                ! Need to specify as much core values are there are
                ! flux surfaces
                if (size(options%ffvoptions%outerfluxval, 1) /= fscc) then
                    print *, 'number of outer flux surfaces: ', fscc-ncfs
                    print *, 'field line IDs of core flux surfaces: ', fsind(ncfs+1:fscc)
                    call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                        'need to specify the amount of outer flux surface values mentioned above')
                end if

                ! Add
                psid(ncfs+1:fscc) = options%ffvoptions%outerfluxval

            case default

                ! Throw error
                call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                    'unknown method to determine flux values for outer flux')

            end select

        end if

        ! Update
        ncfs = fscc

        ! Lonely flux surfaces
        !=====================
        ! Basically flux surfaces with only one non-zero neighbour
        if (dotpflux) then
            ! Loop over all flux surfaces
            do i = 1, nfs
                ! Initialize
                islonely = .true.

                ! Skip if already constrained
                if (hasbeenfound(i)) then
                    islonely = .false.
                    cycle
                end if

                ! Get all vertices of this flux surface
                allocate(tv(count(grid%vert%fieldlineID == i)))
                tv = pack(allvertIDs, grid%vert%fieldlineID == i)

                ! Hedge for no vertices
                if (size(tv) == 0) then
                    deallocate(tv)
                    islonely = .false.
                    cycle
                end if

                ! Loop over all vertices
                do j = 1, size(tv)
                    ! Get neighbours & fieldline IDs of this vertex
                    tvn = GetVertNeig(grid%vert, tv(j))

                    ! If all IDs are either equal to the current flux surface ID
                    ! or equal to zero, or equal to a tangency point ID,
                    ! or equal to the unique flux surface ID that is not
                    ! zero and not the current ID, then continue.
                    ! Otherwise, this flux surface is not lonely
                    tvnID = fieldlineID(tvn)
                    vID = fieldlineID(tv(j))
                    call Unique(pack(tvnID, tvnID /= 0 .and. tvnID /= vID), uniqueID)
                    if (size(uniqueID) == 1) then
                        if (any( (tvnID /= vID) .and. (tvnID /= 0) .and. &
                            (tvnID /= uniqueID(1)) .and. (.not. istangencypoint(tvn)))) then
                            islonely = .false.
                            exit
                        end if
                    else
                        islonely = .false.
                        exit
                    end if
                end do

                ! Add if not lonely
                if (islonely .and. (size(tv) > 2) .and. (.not. any(istangencypoint(tv)))) then
                    ! Add
                    fsind(fscc+1) = i
                    hasbeenfound(i) = .true.

                    ! Update counter
                    fscc = fscc + 1
                end if

                ! Housekeeping
                deallocate(tv)
            end do

            ! Determine flux value
            select case (options%ffvoptions%fixouterfluxmeth)

            case ('auto')

                ! Precompute flux values
                call magneticField%interp%Evaluate(x, y, 0, 0, psi)

                ! Precompute vessel vertices
                call DetermineVesselVertices(isvesselvertex, isvesselface, grid)

                ! Determine flux value as mean of current flux values
                do i = ncfs+1, fscc
                    nfsv = 0
                    do j = 1, grid%vert%ntot
                        if ( (fieldlineID(j) == fsind(i)) .and. (isvesselvertex(j))) then

                            ! Compute
                            psid(i) = psid(i) + psi(j)
                            nfsv = nfsv + 1

                        end if
                    end do

                    if (nfsv == 0) then
                        ! No vertices found, throw error
                        print *, 'outer flux surface ID: ', fsind(i)
                        call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                            'no boundary vertices found for outer flux surface ID')
                    end if

                    ! Average
                    psid(i) = psid(i)/nfsv

                end do

            case ('manual')

                call gdErrorHandler('InitializeFixedFluxValuesConstraints: ' // &
                    'manual option not available for lonely flux surfaces')

            case default

                ! Throw error
                call gdErrorHandler('InitializeFixedFluxvaluesConstraints: ' // &
                    'unknown method to determine flux values for outer flux')

            end select

        end if

        ! Update constraint counter and issue warning if necessary
        dowarning  = .false.
        do i = 1, grid%vert%ntot
            if (any(fieldlineID(i) == fsind(1:fscc))) then
                ccv(i) = ccv(i) + 1
                if (ccv(i) > maxccv(i)) then
                    dowarning = .true.
                end if
            end if
        end do

        if (dowarning) then
            print *, 'InitializeFixedFluxvaluesConstraints: by imposing ' // &
                'fixed psi value, some vertices may be overly constrained!'
        end if
        ! Add
        !====
        constraints%psid = psid(1:fscc)
        constraints%fsind = fsind(1:fscc)
        constraints%ncon = fscc
        allocate(constraints%psiind(fscc))
        constraints%psiind = 0

        ! Write data
        !===========
        if (options%writedata == 1) then
            ! Loop over all constrained flux surfaces and get vertices
            if (allocated(tv)) deallocate(tv)
            allocate(tv(0))
            do i = 1, constraints%ncon
                tv = [tv, pack(allvertIDs, grid%vert%fieldlineID == constraints%fsind(i))]
            end do
            call WriteVertexData(tv, grid%vert%x(tv), grid%vert%y(tv), 'con_ffv_vertices')
        end if

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Evaluation
    subroutine EvaluateFixedFluxvaluesConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! The flux value constraints are evaluated per flux surface
        ! index.


        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(FixedFluxvaluesConstraintsUDT)    :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, k
        integer(I8), allocatable            :: valindex(:), conindex(:), &
            psiind(:)

        ! Auxiliary variables

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if
        if (present(dGdvarin)) then
            dGdvar = dGdvarin
        end if
        if (present(dgradGdvarin)) then
            dgradGdvar = dgradGdvarin
        end if


        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select

        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            interp  => magneticField%interp,    &
            nc      => constraints%ncon,        &
            vert    => grid%vert,               &
            nv      => grid%vert%ntot,          &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Unpack
        select type (designvariables)

        type is (DesignVariablesCoordinatesFluxUDT)

            ! Unpack
            allocate(psiind, source=designvariables%psiind)
            psiind = designvariables%psiind

        class default

            ! Throw error
            call gdErrorHandler('EvaluateFixedFluxvaluesConstraints: ' // &
                'design variable type should be coordinates_desiredflux')

        end select

        ! Constraint value
        !=================
        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateFixedFluxvaluesConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Compute
        G = designvariables%phi(constraints%psiind) - constraints%psid

        ! Constraint gradient
        !====================
        if (dogradient) then
            ! Initialize
            jacG%nrow = nc
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates_desiredflux') ! only flux contributions

                ! Order in jacobian: first x, then y.

                ! Allocate
                jacG%nval = nc
                call jacG%Allocate()
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! psi-contribution
                !---------------
                ! Build indices for xv1
                valindex = [(k, k = ivg+1, ivg+nc)]

                ! Add values
                jacG%row(valindex) = conindex
                jacG%col(valindex) = constraints%psiind
                jacG%val(valindex) = 1

                ! Update counters
                ivg = ivg + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol
                gradG%ncol = jacG%nrow
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        if (dohessian) then

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi
            hessG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates_desiredflux') ! only flux, but no contributions

                ! Initialize
                !===========
                ! Allocate
                hessG%nval = 0
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(nc))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if


            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                      LINE FOLDING CONSTRAINTS                    !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeLinefoldingConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the line folding constraints. There are three
        ! line folding 'types': poloidal, radial, and vessel. Each of
        ! these depends on a different coordinate field: poloidal is
        ! simply the poloidal magnetic field vector, radial a direction
        ! perpendicular in the 2D plane, and vessel is the coordinate
        ! direction that follows vessel contours.

        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(LinefoldingConstraintsUDT)            :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                 :: designvariables

        ! Auxiliary
        integer(I8)                                 :: nvpairspol, &
            nvpairsrad, nvpairsves
        integer(I8), allocatable                    :: &
            vpairspol(:, :), vpairsrad(:, :), vpairsves(:, :), tvID(:), &
            tv(:)

        logical, allocatable                        :: isvesselvertex(:), &
            isvesselface(:)

        real(R8), allocatable                       :: xf(:), yf(:), &
            xv(:, :), yv(:, :), gx(:), gy(:), myones(:), dx(:), dy(:), &
            dotprod(:), signvecpol(:), signvecrad(:), signvecves(:)

        ! Loop
        integer(I8)                                 :: i, j

        ! Initialize
        !===========
        ! Set constants
        constraints%smallnumber = options%lfoptions%smallnumber

        ! Associate
        associate(&
            maxccv          => monitor%maxeqvcc,        &
            ccv             => monitor%eqvcc,           &
            vert            => grid%vert,               &
            face            => grid%face,               &
            x               => grid%vert%x,             &
            y               => grid%vert%y,             &
            fieldlineID     => grid%vert%fieldlineID    &
            )

        ! Determine vessel faces
        call DetermineVesselVertices(isvesselvertex, isvesselface, grid)

        ! Determine poloidal pairs
        !=========================
        nvpairspol = 0
        allocate(vpairspol(face%ntot, 2))
        vpairspol = 0
        if (options%lfoptions%poloidal) then
            do i = 1, face%ntot
                ! Skip vessel faces
                if (isvesselface(i)) then
                    cycle
                end if

                ! Get vertices
                tv = face%vert(i, :)

                ! Get IDs
                tvID = fieldlineID(tv)

                ! If zero or not the same, continue
                if (any(tvID == 0)) then
                    cycle
                end if
                if (tvID(1) /= tvID(2)) then
                    cycle
                end if

                ! Add vertex pair
                nvpairspol = nvpairspol + 1
                vpairspol(nvpairspol, :) = tv(1:2)

            end do
        end if

        ! Trim
        vpairspol = vpairspol(1:nvpairspol, :)

        ! Determine radial pairs
        !=======================
        nvpairsrad = 0
        allocate(vpairsrad(face%ntot, 2))
        vpairsrad = 0
        if (options%lfoptions%radial) then
            do i = 1, face%ntot
                ! Skip vessel faces
                if (isvesselface(i)) then
                    cycle
                end if

                ! Get vertices
                tv = face%vert(i, :)

                ! Get IDs
                tvID = fieldlineID(tv)

                ! If zero or the same, continue
                if (any(tvID == 0)) then
                    cycle
                end if
                if (tvID(1) == tvID(2)) then
                    cycle
                end if

                ! Add vertex pair
                nvpairsrad = nvpairsrad + 1
                vpairsrad(nvpairsrad, :) = tv(1:2)
            end do
        end if

        ! Trim
        vpairsrad = vpairsrad(1:nvpairsrad, :)

        ! Determine vessel pairs
        !=======================
        nvpairsves = 0
        allocate(vpairsves(face%ntot, 2))
        vpairsves = 0
        if (options%lfoptions%vessel) then
            do i = 1, face%ntot
                ! Take only vessel faces
                if (.not. isvesselface(i)) then
                    cycle
                end if

                ! Get vertices
                tv = face%vert(i, :)

                ! Get IDs
                tvID = fieldlineID(tv)

                ! At least one vertex should have zero ID
                if (.not. any(tvID == 0)) then
                    cycle
                end if

                ! Add vertex pair
                nvpairsves = nvpairsves + 1
                vpairsves(nvpairsves, :) = tv(1:2)
            end do
        end if

        ! Trim
        vpairsves = vpairsves(1:nvpairsves, :)

        ! Check orientation
        !==================
        ! Poloidal
        !---------
        ! Allocate
        allocate(xv(nvpairspol, 2), yv(nvpairspol, 2), gx(nvpairspol), &
            gy(nvpairspol), myones(nvpairspol))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairspol(:, j))
            yv(:, j) = y(vpairspol(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call magneticField%interp%Evaluate(xf, yf, 1, 0, gx)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = -gy*dx + gx*dy
        signvecpol = sign(myones, dotprod)

        ! Write
        call WriteVertexPairData(vpairspol, xv, yv, 'con_lf_vpairspol')

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Radial
        !-------
        ! Allocate
        allocate(xv(nvpairsrad, 2), yv(nvpairsrad, 2), gx(nvpairsrad), &
            gy(nvpairsrad), myones(nvpairsrad))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairsrad(:, j))
            yv(:, j) = y(vpairsrad(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call magneticField%interp%Evaluate(xf, yf, 1, 0, gx)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = gx*dx + gy*dy
        signvecrad = sign(myones, dotprod)

        ! Write
        call WriteVertexPairData(vpairsrad, xv, yv, 'con_lf_vpairsrad')

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Vessel
        !-------
        ! Allocate
        allocate(xv(nvpairsves, 2), yv(nvpairsves, 2), gx(nvpairsves), &
            gy(nvpairsves), myones(nvpairsves))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairsves(:, j))
            yv(:, j) = y(vpairsves(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call environment%vessel%plfvessel%Evaluate(xf, yf, 1, 0, gx)
        call environment%vessel%plfvessel%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = -gy*dx + gx*dy
        signvecves = sign(myones, dotprod)

        ! Write
        call WriteVertexPairData(vpairsves, xv, yv, 'con_lf_vpairsves')

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Add
        !----
        constraints%vpairspol   = vpairspol
        constraints%vpairsrad   = vpairsrad
        constraints%vpairsves   = vpairsves

        constraints%nvpairspol  = nvpairspol
        constraints%nvpairsrad  = nvpairsrad
        constraints%nvpairsves  = nvpairsves

        constraints%signvecpol  = signvecpol
        constraints%signvecrad  = signvecrad
        constraints%signvecves  = signvecves

        constraints%ncon        = nvpairspol + nvpairsrad + nvpairsves


        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! Evaluation
    subroutine EvaluateLinefoldingConstraints(constraints, G, gradG, &
        hessG, grid, magneticField, environment, dogradient, &
        dohessian, designvariables, lambda, varin, valuesin, dGdvarin, &
        dgradGdvarin)

        ! Description
        !============
        ! Evaluate the line folding constraints. These constraints are formulated
        ! as follows per face considered:
        !
        !           H(i) = -sign(dotprod(B, v)*signvec + tol <= 0.
        !
        ! Here, B is the (magnetic) field that defines the coordinate line
        ! direction, 'signvec' is a predetermined multiplier, v is the vector
        ! from the first point to the second of the vertex pairs considered, vx and
        ! vy are the x- and y-components of this vector, and tol is a constant and
        ! small number that prevents vertex coincidence (and, if alignement
        ! constraints are imposed, a measure of the minimal distance along the
        ! coordinate line). This formulation of the constraint has the following
        ! advantages:
        ! - the gradient w.r.t. the coordinates exists always, even if vertices
        ! coincide (vx = vy = 0).
        ! - vertices are projected along the correct direction when the constraint
        ! is violated
        !
        ! The following downsides require caution:
        ! - if the (magnetic) field exhibits nulls, the constraint cannot be
        ! properly imposed. This is checked and hedged for (i.e. the constraint is
        ! simply set to be inactive...). Warnings will be issued.
        ! - We need to hedge to prevent imposing the same constraint multiple
        ! times (see below) to prevent infeasible subproblems
        ! - We need to assume that the initial signvec vector is properly chosen,
        ! otherwise this constraint will enforce linefolding. This typically
        ! requires proper initial knowledge of the grid.
        ! - The tolerance should be chosen appropriately to the grid dimensions and
        ! should not interfere with other constraints.

        ! Notes
        !======

        ! Declare variables
        !==================
        ! Arguments
        class(LinefoldingConstraintsUDT)    :: constraints
        real(R8), allocatable               :: G(:)
        real(R8), allocatable               :: lambda(:)
        type(MySparseUDT)                   :: hessG, gradG, jacG
        type(GridUDT)                       :: grid
        type(MagneticFieldUDT)              :: magneticField
        type(EnvironmentUDT)                :: environment
        logical                             :: dogradient, dohessian
        class(DesignVariablesGDUDT)         :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8), intent(in), optional      :: valuesin(:)
        type(MySparseUDT), optional         :: dGdvarin, dgradGdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        type(MySparseUDT)                   :: dGdvar, dgradGdvar

        ! Auxiliary
        integer(I8)                             :: nvpairs, ntvp
        integer(I8), allocatable, dimension(:)  :: tvp

        real(R8), allocatable, dimension(:, :)  :: xvp, yvp, xvr, yvr, &
            xvv, yvv
        real(R8), allocatable, dimension(:)     :: xfp, yfp, dxp, dyp, &
            gxp, gyp, xfr, yfr, dxr, dyr, gxr, gyr, xfv, yfv, dxv, dyv, &
            gxv, gyv, gxxp, gxxr, gxxv, gxyp, gxyr, gxyv, gyxp, gyxr, &
            gyxv, gyyp, gyyr, gyyv, gxx, gxy, gyx, gyy, gxxx, gxxy, &
            gxyx, gyxx, gxyy, gyxy, gyyx, gyyy, gxxxp, gxxyp, &
            gxyxp, gyxxp, gxyyp, gyxyp, gyyxp, gyyyp, gxxxr, gxxyr, &
            gxyxr, gyxxr, gxyyr, gyxyr, gyyxr, gyyyr, gxxxv, gxxyv, &
            gxyxv, gyxxv, gxyyv, gyxyv, gyyxv, gyyyv, xf, yf, dx, dy, &
            gx, gy, signvec

        ! Loop variables
        integer(I8)                         :: ic, ivg, ivh, j, k
        integer(I8), allocatable            :: valindex(:), conindex(:), &
            vpairs(:, :)

        ! Initialize
        !===========
        ! Check inputs
        if (present(varin)) then
            var = varin
        else
            var = 'no'
        end if
        if (present(valuesin)) then
            values = valuesin
        else
            allocate(values(0))
        end if

        ! Checks
        if ( (.not. allocated(lambda)) .and. dohessian) then
            ! Throw error
            call gdErrorHandler('When evaluating the hessian vector' &
                // ' multiplication, lambda must be given')
        end if

        if (size(lambda) .ne. constraints%ncon) then
            ! Lambda should have the same size as the constraints
            call gdErrorHandler('Lambda should have the same size ' &
                // 'as the constraint vector')
        end if

        ! Check derivative computation
        select case (var)

        case ('no')

            ! No derivatives, initialize correctly
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case ('vesselcoordinates')

            ! No contributions
            dGdvar = SpZeros(constraints%ncon, size(values, 1))
            dgradGdvar = SpZeros(designvariables%nphi, size(values, 1))

        case default

            ! Not implemented
            call gdErrorHandler('EvaluateFluxfunctionConstraints: variable not implemented')

        end select
        ! Counters
        ic = 0 ! constraint counter (local)
        ivg = 0 ! value index for gradient
        ivh = 0 ! value index for hessian

        ! Associate
        associate(&
            signvecpol      => constraints%signvecpol,  &
            signvecrad      => constraints%signvecrad,  &
            signvecves      => constraints%signvecves,  &
            vpairspol       => constraints%vpairspol,   &
            vpairsrad       => constraints%vpairsrad,   &
            vpairsves       => constraints%vpairsves,   &
            nvpairspol      => constraints%nvpairspol,  &
            nvpairsrad      => constraints%nvpairsrad,  &
            nvpairsves      => constraints%nvpairsves,  &
            tol             => constraints%smallnumber, &
            nc      => constraints%ncon,        &
            x       => grid%vert%x,             &
            y       => grid%vert%y              &
            )

        ! Allocate
        if (.not. allocated(G)) then
            allocate(G(nc))
        else
            if (size(G) .ne. nc) then

                ! Print a warning and reallocate
                print *, 'EvaluateLinefoldingConstraints: ' &
                    // 'Wrong dimension of G, reallocating'

                ! Deallocate and reallocate
                deallocate(G)
                allocate(G(nc))

            end if
        end if

        ! Constraint value
        !=================
        ! Precompute
        !-----------
        ! Allocate
        allocate(xvp(nvpairspol, 2), yvp(nvpairspol, 2), &
            xfp(nvpairspol), yfp(nvpairspol), dxp(nvpairspol), &
            dyp(nvpairspol), gxp(nvpairspol), gyp(nvpairspol))
        allocate(xvr(nvpairsrad, 2), yvr(nvpairsrad, 2), &
            xfr(nvpairsrad), yfr(nvpairsrad), dxr(nvpairsrad), &
            dyr(nvpairsrad), gxr(nvpairsrad), gyr(nvpairsrad))
        allocate(xvv(nvpairsves, 2), yvv(nvpairsves, 2), &
            xfv(nvpairsves), yfv(nvpairsves), dxv(nvpairsves), &
            dyv(nvpairsves), gxv(nvpairsves), gyv(nvpairsves))

        ! Poloidal
        !---------
        ! Precompute
        do j = 1, 2
            xvp(:, j) = x(vpairspol(:, j))
            yvp(:, j) = y(vpairspol(:, j))
        end do
        dxp = xvp(:, 2) - xvp(:, 1)
        dyp = yvp(:, 2) - yvp(:, 1)
        xfp = 0.5*(xvp(:, 1) + xvp(:, 2))
        yfp = 0.5*(yvp(:, 1) + yvp(:, 2))
        call magneticField%interp%Evaluate(xfp, yfp, 1, 0, gyp)
        call magneticField%interp%Evaluate(xfp, yfp, 0, 1, gxp)
        gxp = -gxp

        ! Evaluate
        G(ic+1:ic+nvpairspol) = -signvecpol*(gxp*dxp + gyp*dyp) + tol
        if (any(G(ic+1:ic+nvpairspol) > 0)) then
            !print *, 'EvaluateLineFoldingConstraints: poloidal folding active'
        end if
        ntvp = count(G(ic+1:ic+nvpairspol) > 0)
        allocate(tvp(ntvp))
        tvp = pack([(k, k = 1, nvpairspol)], G(ic+1:ic+nvpairspol) > 0)
        call WriteVertexPairData(vpairspol(tvp, :), xvp(tvp, :), yvp(tvp, :), &
            'con_lf_vpairspol_iterate')
        deallocate(tvp)

        ! Update counter
        ic = ic + nvpairspol

        ! Radial
        !-------
        ! Precompute
        do j = 1, 2
            xvr(:, j) = x(vpairsrad(:, j))
            yvr(:, j) = y(vpairsrad(:, j))
        end do
        dxr = xvr(:, 2) - xvr(:, 1)
        dyr = yvr(:, 2) - yvr(:, 1)
        xfr = 0.5*(xvr(:, 1) + xvr(:, 2))
        yfr = 0.5*(yvr(:, 1) + yvr(:, 2))
        call magneticField%interp%Evaluate(xfr, yfr, 1, 0, gxr)
        call magneticField%interp%Evaluate(xfr, yfr, 0, 1, gyr)

        ! Evaluate
        G(ic+1:ic+nvpairsrad) = -signvecrad*(gxr*dxr + gyr*dyr) + tol
        if (any(G(ic+1:ic+nvpairsrad) > 0)) then
            !print *, 'EvaluateLineFoldingConstraints: radial folding active'
        end if
        ntvp = count(G(ic+1:ic+nvpairsrad) > 0)
        allocate(tvp(ntvp))
        tvp = pack([(k, k = 1, nvpairsrad)], G(ic+1:ic+nvpairsrad) > 0)
        call WriteVertexPairData(vpairsrad(tvp, :), xvr(tvp, :), yvr(tvp, :), &
            'con_lf_vpairsrad_iterate')
        deallocate(tvp)

        ! Update counter
        ic = ic + nvpairsrad

        ! Vessel
        !-------
        ! Precompute
        do j = 1, 2
            xvv(:, j) = x(vpairsves(:, j))
            yvv(:, j) = y(vpairsves(:, j))
        end do
        dxv = xvv(:, 2) - xvv(:, 1)
        dyv = yvv(:, 2) - yvv(:, 1)
        xfv = 0.5*(xvv(:, 1) + xvv(:, 2))
        yfv = 0.5*(yvv(:, 1) + yvv(:, 2))
        call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 0, gyv)
        call environment%vessel%plfvessel%Evaluate(xfv, yfv, 0, 1, gxv)
        gxv = -gxv

        ! Evaluate
        G(ic+1:ic+nvpairsves) = -signvecves*(gxv*dxv + gyv*dyv) + tol
        if (any(G(ic+1:ic+nvpairsves) > 0)) then
            !print *, 'EvaluateLineFoldingConstraints: vessel folding active'
        end if
        ntvp = count(G(ic+1:ic+nvpairsves) > 0)
        allocate(tvp(ntvp))
        tvp = pack([(k, k = 1, nvpairsves)], G(ic+1:ic+nvpairsves) > 0)
        call WriteVertexPairData(vpairsves(tvp, :), xvv(tvp, :), yvv(tvp, :), &
            'con_lf_vpairsves_iterate')
        deallocate(tvp)


        ! Update counter
        ic = ic + nvpairsves

        ! Constraint gradient
        !====================
        ! Reset counter
        ic = 0

        ! Precompute
        if (dohessian .or. dogradient) then

            ! Poloidal
            !---------
            ! Allocate
            allocate(gxxp(nvpairspol), gxyp(nvpairspol), &
            gyxp(nvpairspol), gyyp(nvpairspol))

            ! Precompute
            call magneticField%interp%Evaluate(xfp, yfp, 1, 1, gxxp)
            call magneticField%interp%Evaluate(xfp, yfp, 0, 2, gxyp)
            call magneticField%interp%Evaluate(xfp, yfp, 2, 0, gyxp)
            call magneticField%interp%Evaluate(xfp, yfp, 1, 1, gyyp)
            gxxp = -gxxp
            gxyp = -gxyp

            ! Radial
            !-------
            ! Allocate
            allocate(gxxr(nvpairsrad), gxyr(nvpairsrad), &
            gyxr(nvpairsrad), gyyr(nvpairsrad))

            ! Precompute
            call magneticField%interp%Evaluate(xfr, yfr, 2, 0, gxxr)
            call magneticField%interp%Evaluate(xfr, yfr, 1, 1, gxyr)
            call magneticField%interp%Evaluate(xfr, yfr, 1, 1, gyxr)
            call magneticField%interp%Evaluate(xfr, yfr, 0, 2, gyyr)

            ! Vessel
            !-------
            ! Allocate
            allocate(gxxv(nvpairsves), gxyv(nvpairsves), &
            gyxv(nvpairsves), gyyv(nvpairsves))

            ! Precompute
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 1, gxxv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 0, 2, gxyv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 2, 0, gyxv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 1, gyyv)
            gxxv = -gxxv
            gxyv = -gxyv

            ! Concatenate
            !------------
            gxx = [gxxp, gxxr, gxxv]
            gxy = [gxyp, gxyr, gxyv]
            gyx = [gyxp, gyxr, gyxv]
            gyy = [gyyp, gyyr, gyyv]
            gx = [gxp, gxr, gxv]
            gy = [gyp, gyr, gyv]

            xf = [xfp, xfr, xfv]
            yf = [yfp, yfr, yfv]
            dx = [dxp, dxr, dxv]
            dy = [dyp, dyr, dyv]

            nvpairs = nvpairspol + nvpairsrad + nvpairsves
            allocate(vpairs(nvpairs, 2))
            do j = 1, 2
                vpairs(:, j) = [vpairspol(:, j), vpairsrad(:, j), vpairsves(:, j)]
            end do
            signvec = [signvecpol, signvecrad, signvecves]

        end if

        if (dogradient) then
            ! Initialize
            jacG%nrow = nc
            jacG%ncol = designvariables%nphi

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux') ! no flux contributions

                ! Order in jacobian: first x, then y.

                ! Allocate
                jacG%nval = 4*nc ! 4 contributions per constraint
                call jacG%Allocate()
                allocate(conindex(nc))
                allocate(valindex(nc))

                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! Add values
                valindex = [(k, k = ivg+1, ivg+nc)]
                jacG%row(valindex) = conindex
                jacG%col(valindex) = vpairs(:, 1)
                jacG%val(valindex) = -signvec*(-gx + 0.5*dx*gxx + 0.5*dy*gyx) ! x1
                ivg = ivg + nc

                valindex = valindex + nc
                jacG%row(valindex) = conindex
                jacG%col(valindex) = vpairs(:, 1) + grid%vert%ntot
                jacG%val(valindex) = -signvec*(-gy + 0.5*dx*gxy + 0.5*dy*gyy) !y1
                ivg = ivg + nc

                valindex = valindex + nc
                jacG%row(valindex) = conindex
                jacG%col(valindex) = vpairs(:, 2)
                jacG%val(valindex) = -signvec*(gx + 0.5*dx*gxx + 0.5*dy*gyx) ! x2
                ivg = ivg + nc

                valindex = valindex + nc
                jacG%row(valindex) = conindex
                jacG%col(valindex) = vpairs(:, 2) + grid%vert%ntot
                jacG%val(valindex) = -signvec*(gy + 0.5*dy*gyy + 0.5*dx*gxy) !y2
                ivg = ivg + nc

                ! Update
                ic = ic + nc

                ! Build gradient
                !---------------
                gradG%nrow = jacG%ncol
                gradG%ncol = jacG%nrow
                gradG%nval = jacG%nval

                call gradG%Allocate()
                gradG%row = jacG%col
                gradG%col = jacG%row
                gradG%val = jacG%val

                ! Housekeeping
                call jacG%Deallocate()

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Constraint hessian
        !===================
        ! Reset counter
        ic = 0
        if (dohessian) then

            ! Initialize
            ic = 0
            hessG%nrow = designvariables%nphi
            hessG%ncol = designvariables%nphi

            ! Precompute
            !===========
            ! Poloidal
            !---------
            ! Allocate
            allocate(gxxxp(nvpairspol), gxyxp(nvpairspol), &
                gyxxp(nvpairspol), gyyxp(nvpairspol), gxxyp(nvpairspol), &
                gxyyp(nvpairspol), gyxyp(nvpairspol), gyyyp(nvpairspol))

            ! Precompute
            call magneticField%interp%Evaluate(xfp, yfp, 2, 1, gxxxp)
            call magneticField%interp%Evaluate(xfp, yfp, 1, 2, gxyxp)
            call magneticField%interp%Evaluate(xfp, yfp, 3, 0, gyxxp)
            call magneticField%interp%Evaluate(xfp, yfp, 2, 1, gyyxp)

            call magneticField%interp%Evaluate(xfp, yfp, 1, 2, gxxyp)
            call magneticField%interp%Evaluate(xfp, yfp, 0, 3, gxyyp)
            call magneticField%interp%Evaluate(xfp, yfp, 2, 1, gyxyp)
            call magneticField%interp%Evaluate(xfp, yfp, 1, 2, gyyyp)

            gxxxp = -gxxxp
            gxyxp = -gxyxp
            gxxyp = -gxxyp
            gxyyp = -gxyyp

            ! Radial
            !-------
            ! Allocate
            allocate(gxxxr(nvpairsrad), gxyxr(nvpairsrad), &
                gyxxr(nvpairsrad), gyyxr(nvpairsrad), gxxyr(nvpairsrad), &
                gxyyr(nvpairsrad), gyxyr(nvpairsrad), gyyyr(nvpairsrad))

            ! Precompute
            call magneticField%interp%Evaluate(xfr, yfr, 3, 0, gxxxr)
            call magneticField%interp%Evaluate(xfr, yfr, 2, 1, gxyxr)
            call magneticField%interp%Evaluate(xfr, yfr, 2, 1, gyxxr)
            call magneticField%interp%Evaluate(xfr, yfr, 1, 2, gyyxr)

            call magneticField%interp%Evaluate(xfr, yfr, 2, 1, gxxyr)
            call magneticField%interp%Evaluate(xfr, yfr, 1, 2, gxyyr)
            call magneticField%interp%Evaluate(xfr, yfr, 1, 2, gyxyr)
            call magneticField%interp%Evaluate(xfr, yfr, 0, 3, gyyyr)

            ! Vessel
            !-------
            ! Allocate
            allocate(gxxxv(nvpairsves), gxyxv(nvpairsves), &
                gyxxv(nvpairsves), gyyxv(nvpairsves), gxxyv(nvpairsves), &
                gxyyv(nvpairsves), gyxyv(nvpairsves), gyyyv(nvpairsves))

            ! Precompute
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 2, 1, gxxxv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 2, gxyxv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 3, 0, gyxxv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 2, 1, gyyxv)

            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 2, gxxyv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 0, 3, gxyyv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 2, 1, gyxyv)
            call environment%vessel%plfvessel%Evaluate(xfv, yfv, 1, 2, gyyyv)

            gxxxv = -gxxxv
            gxyxv = -gxyxv
            gxxyv = -gxxyv
            gxyyv = -gxyyv

            ! Concatenate
            !------------
            gxxx = [gxxxp, gxxxr, gxxxv]
            gxyx = [gxyxp, gxyxr, gxyxv]
            gyxx = [gyxxp, gyxxr, gyxxv]
            gyyx = [gyyxp, gyyxr, gyyxv]

            gxxy = [gxxyp, gxxyr, gxxyv]
            gxyy = [gxyyp, gxyyr, gxyyv]
            gyxy = [gyxyp, gyxyr, gyxyv]
            gyyy = [gyyyp, gyyyr, gyyyv]

            ! Check design variables
            select case(designvariables%type)

            case ('coordinates', 'coordinates_desiredflux') ! no flux, only coordinates

                ! Initialize
                !===========
                ! Allocate
                hessG%nval = 16*nc
                if (.not. allocated(valindex)) then
                    allocate(valindex(nc))
                end if
                if (.not. allocated(conindex)) then
                    allocate(conindex(nc))
                end if
                if (.not. allocated(hessG%val)) then
                    call hessG%Allocate()
                end if

                ! Concatenate to reduce lines of code...


                ! Poloidal
                !---------
                ! Build constraint indices
                conindex = [(k, k = ic+1, ic+nc)]

                ! Add values
                ! x1x1
                valindex = [(k, k = ivh+1, ivh+nvpairs)]
                hessG%row(valindex) = vpairs(:, 1)
                hessG%col(valindex) = vpairs(:, 1)
                hessG%val(valindex) = -signvec*(0.25*dx*gxxx - 1.0*gxx + 0.25*dy*gyxx)*lambda(conindex) ! x1x1
                ivh = ivh + nvpairs

                ! x1y1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1)
                hessG%col(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.25*dx*gxyx - 0.5*gyx - 0.5*gxy + 0.25*dy*gyyx)*lambda(conindex) ! x1y1
                ivh = ivh + nvpairs

                ! y1x1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 1)
                hessG%val(valindex) = -signvec*(0.25*dx*gxyx - 0.5*gyx - 0.5*gxy + 0.25*dy*gyyx)*lambda(conindex) ! x1y1
                ivh = ivh + nvpairs

                ! y1y1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.25*dx*gxyy - 1.0*gyy + 0.25*dy*gyyy)*lambda(conindex) ! y1y1
                ivh = ivh + nvpairs

                ! x1x2
                valindex = [(k, k = ivh+1, ivh+nvpairs)]
                hessG%row(valindex) = vpairs(:, 1)
                hessG%col(valindex) = vpairs(:, 2)
                hessG%val(valindex) = -signvec*(0.25*dx*gxxx + 0.25*dy*gyxx)*lambda(conindex) ! x1x2
                ivh = ivh + nvpairs

                ! x1y2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1)
                hessG%col(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.5*gyx - 0.5*gxy + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! x1y2
                ivh = ivh + nvpairs

                ! y1x2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 2)
                hessG%val(valindex) = -signvec*(0.5*gxy - 0.5*gyx + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! y1x2
                ivh = ivh + nvpairs

                ! y1y2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.25*dx*gxyy + 0.25*dy*gyyy)*lambda(conindex) ! y1y2
                ivh = ivh + nvpairs

                ! x2x1
                valindex = [(k, k = ivh+1, ivh+nvpairs)]
                hessG%row(valindex) = vpairs(:, 2)
                hessG%col(valindex) = vpairs(:, 1)
                hessG%val(valindex) = -signvec*(0.25*dx*gxxx + 0.25*dy*gyxx)*lambda(conindex) ! x1x2
                ivh = ivh + nvpairs

                ! x2y1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2)
                hessG%col(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.5*gxy - 0.5*gyx + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! y1x2
                ivh = ivh + nvpairs

                ! y2x1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 1)
                hessG%val(valindex) = -signvec*(0.5*gyx - 0.5*gxy + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! x1y2
                ivh = ivh + nvpairs

                ! y2y1
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 1) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.25*dx*gxyy + 0.25*dy*gyyy)*lambda(conindex) ! y1y2
                ivh = ivh + nvpairs

                ! x2x2
                valindex = [(k, k = ivh+1, ivh+nvpairs)]
                hessG%row(valindex) = vpairs(:, 2)
                hessG%col(valindex) = vpairs(:, 2)
                hessG%val(valindex) = -signvec*(1.0*gxx + 0.25*dx*gxxx + 0.25*dy*gyxx)*lambda(conindex) ! x2x2
                ivh = ivh + nvpairs

                ! x2y2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2)
                hessG%col(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(0.5*gxy + 0.5*gyx + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! x2y2
                ivh = ivh + nvpairs

                ! y2x2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 2)
                hessG%val(valindex) = -signvec*(0.5*gxy + 0.5*gyx + 0.25*dx*gxyx + 0.25*dy*gyyx)*lambda(conindex) ! x2y2
                ivh = ivh + nvpairs

                ! y2y2
                valindex = valindex + nvpairs
                hessG%row(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%col(valindex) = vpairs(:, 2) + grid%vert%ntot
                hessG%val(valindex) = -signvec*(1.0*gyy + 0.25*dx*gxyy + 0.25*dy*gyyy)*lambda(conindex) ! y2y2
                ivh = ivh + nvpairs

                ! Update
                ic = ic + nvpairs

            case default

                ! Unknown, throw error
                call gdErrorHandler('Gradient not implemented for ' &
                    // 'this type of design variable')

            end select

        end if

        ! Housekeeping
        !=============
        ! End associate
        end associate

        ! Optional arguments
        if (present(dGdvarin)) then
            dGdvarin = dGdvar
        end if
        if (present(dgradGdvarin)) then
            dgradGdvarin = dgradGdvar
        end if

    end subroutine

    ! Update
    subroutine UpdateLinefoldingConstraints(constraints, grid, &
        magneticField, environment)

        ! Description
        !============
        ! Update the linefolding constraints description according to the
        ! given grid, magnetic field and environment. Can be used to
        ! update constraint parameters after external updating of these
        ! quantities (e.g. when doing shape optimization, vessel will
        ! change etc).

        ! Note: nothing has to be updated here (yet), since environment
        ! etc should be updated elsewhere and are parsed to the
        ! constraint directly

        ! Declare variables
        !==================
        ! Arguments
        class(LinefoldingConstraintsUDT)            :: constraints
        type(GridUDT), intent(in)                   :: grid
        type(MagneticFieldUDT), intent(in)          :: magneticField
        type(EnvironmentUDT), intent(in)            :: environment

        ! Auxiliary
        real(R8), allocatable                       :: xf(:), yf(:), &
            xv(:, :), yv(:, :), gx(:), gy(:), myones(:), dx(:), dy(:), &
            dotprod(:), signvecpol(:), signvecrad(:), signvecves(:)

        ! Loop
        integer(I8)                                 :: j

        ! Initialize
        !===========
        ! Associate
        associate(&
            nvpairspol      => constraints%nvpairspol,  &
            nvpairsrad      => constraints%nvpairsrad,  &
            nvpairsves      => constraints%nvpairsves,  &
            vpairspol       => constraints%vpairspol,   &
            vpairsrad       => constraints%vpairsrad,   &
            vpairsves       => constraints%vpairsves,   &
            vert            => grid%vert,               &
            face            => grid%face,               &
            x               => grid%vert%x,             &
            y               => grid%vert%y,             &
            fieldlineID     => grid%vert%fieldlineID    &
            )

        ! Check orientation
        !==================
        ! Poloidal
        !---------
        ! Allocate
        allocate(xv(nvpairspol, 2), yv(nvpairspol, 2), gx(nvpairspol), &
            gy(nvpairspol), myones(nvpairspol))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairspol(:, j))
            yv(:, j) = y(vpairspol(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call magneticField%interp%Evaluate(xf, yf, 1, 0, gx)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = -gy*dx + gx*dy
        signvecpol = sign(myones, dotprod)

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Radial
        !-------
        ! Allocate
        allocate(xv(nvpairsrad, 2), yv(nvpairsrad, 2), gx(nvpairsrad), &
            gy(nvpairsrad), myones(nvpairsrad))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairsrad(:, j))
            yv(:, j) = y(vpairsrad(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call magneticField%interp%Evaluate(xf, yf, 1, 0, gx)
        call magneticField%interp%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = gx*dx + gy*dy
        signvecrad = sign(myones, dotprod)

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Vessel
        !-------
        ! Allocate
        allocate(xv(nvpairsves, 2), yv(nvpairsves, 2), gx(nvpairsves), &
            gy(nvpairsves), myones(nvpairsves))
        myones = 1

        ! Compute coordinates and vectors
        do j = 1, 2
            xv(:, j) = x(vpairsves(:, j))
            yv(:, j) = y(vpairsves(:, j))
        end do
        dx = xv(:, 2) - xv(:, 1)
        dy = yv(:, 2) - yv(:, 1)
        xf = 0.5*(xv(:, 1) + xv(:, 2))
        yf = 0.5*(yv(:, 1) + yv(:, 2))
        call environment%vessel%plfvessel%Evaluate(xf, yf, 1, 0, gx)
        call environment%vessel%plfvessel%Evaluate(xf, yf, 0, 1, gy)

        ! Evaluate dot product and save sign
        dotprod = -gy*dx + gx*dy
        signvecves = sign(myones, dotprod)

        ! Housekeeping
        deallocate(xv, yv, gx, gy, myones)

        ! Add
        !----
        constraints%signvecpol  = signvecpol
        constraints%signvecrad  = signvecrad
        constraints%signvecves  = signvecves

        ! Housekeeping
        !=============
        end associate


    end subroutine

    !------------------------------------------------------------------!
    !                       IN-VESSEL CONSTRAINTS                      !
    !------------------------------------------------------------------!

    ! Initialize
    subroutine InitializeInVesselConstraints(constraints, &
        grid, magneticField, environment, monitor, designvariables, &
        options)

        ! Description
        !============
        ! Initialize the line folding constraints. There are three
        ! line folding 'types': poloidal, radial, and vessel. Each of
        ! these depends on a different coordinate field: poloidal is
        ! simply the poloidal magnetic field vector, radial a direction
        ! perpendicular in the 2D plane, and vessel is the coordinate
        ! direction that follows vessel contours.

        ! Notes
        !======

        ! Modules
        use mod_definitions, only : targetID, vesselID

        ! Declare variables
        !==================
        ! Arguments
        class(InVesselConstraintsUDT)               :: constraints
        type(GridUDT)                               :: grid
        type(MagneticFieldUDT)                      :: magneticField
        type(EnvironmentUDT)                        :: environment
        type(ConstraintsMonitorUDT)                 :: monitor
        type(ConstraintOptionsUDT)                  :: options
        class(DesignVariablesGDUDT)                 :: designvariables

        ! Auxiliary
        logical, allocatable                        :: isconstrained(:)

        ! Loop
        integer(I8)                                 :: i

        ! Associate
        associate(&
            vessel      => environment%vessel  &
            )

        ! Bookkeeping of constrained vertices (to prevent imposing
        ! constraint twice)
        allocate(isconstrained(grid%vert%ntot))
        isconstrained = .true. ! constrain all, but ignore boundary

        ! Construct boundary
        !===================
        ! Should already be constructed in vessel - assign
        constraints%plf = vessel%plfvessel

        ! Visualize
        call constraints%plf%Visualize('constraints_invessel_plf')

        ! Set the constraints
        !====================
        ! Kick out vertices that are fixed on the boundary
        constraints%ncon = 0
        constraints%nvert = 0
        do i = 1, size(grid%bnd)

            ! Check if target plate - hard coded here...
            if (any(grid%bnd(i)%ID == [targetID, vesselID])) then

                ! Get the current vertices, set to false
                isconstrained(grid%bnd(i)%vert) = .false.

            end if

        end do

        ! Add
        constraints%vert = pack([(i, i = 1, grid%vert%ntot)], isconstrained)
        constraints%nvert = size(constraints%vert)
        constraints%ncon = constraints%nvert

        ! End associate
        end associate

    end subroutine

    ! Update




end module
