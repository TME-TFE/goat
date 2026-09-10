!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the cost function implementation. Each cost
! function inherits from the basic cost function type for the grid 
! deformation, which itself inherits from the archetypical cost function 
! type defined in the optmod_costfunction module. The cost functions 
! defined here are specific for shape optimization. 

module somod_costfunction
    
    ! Initialize
    !============
    ! Load modules
    use optmod_costfunction
    use gdmod_optimizationengine
    use somod_userinput
    use somod_designvariables
    use PolygonLevelsetFunction2D
    use mod_linearsolverinterface
    use optmod_hessianapproximation

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Abstract types
    !===============
    ! General cost function type
    type, abstract, extends(CostfunctionUDT) :: CostfunctionSOUDT

        ! Description
        !============
        ! Defines the basic cost function structure for the shape 
        ! optimization. The following general fields are added: 
        ! - J:          The cost function value (scalar)
        ! - type:       the cost function type (string)
        ! - B:          hessian estimator
        ! - doremesh:   logical flag to indicate if remeshing is 
        !               requested. Cost function and gradient value 
        !               may be garbage in that case. 

        ! The following routines should be implemented for these cost
        ! functions (see also the interface below for a description of
        ! what the routines should do):
        ! - Initialize
        ! - Evaluate

        ! Note that it is up to the developer to use/update the hessian
        ! estimator if required/desired! 

        ! Cost function value
        real(R8)                        :: J 

        ! Cost function type
        character(:), allocatable       :: type

        ! Hessian estimator
        class(HessianApproximationUDT), allocatable     :: B

        ! Remeshing 
        logical                         :: doremesh 

    contains

        ! Cost function initialization
        procedure(InitializeCostfunctionSOINT), deferred      :: Initialize

        ! Cost function evaluation
        procedure(EvaluateCostFunctionSOINT), deferred        :: Evaluate

    end type

    ! Derived types
    !==============
    ! Dummy (zero) cost function (e.g. for reduced formulation with solps)
    type, extends(CostFunctionSOUDT) :: CostFunctionDummyUDT 

        ! Description
        !============
        ! Dummy cost function that simply returns zero  for evaluation
        ! etc.

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeCostFunctionDummy 

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionDummy

    end type
    
    ! Levelset fitting cost function
    type, extends(CostfunctionSOUDT) :: CostfunctionPLFUDT

        ! Description
        !============
        ! This cost function aims to match the vessel coordinates to 
        ! a certain given levelset function. This levelset function 
        ! should be chosen wisely of course. Currently, the 
        ! initialization of this levelset function is based on an 
        ! additional vessel structure file that is read in and for which
        ! a polygon description is made. The following fields are 
        ! present:
        !
        !   targetplf       : target polygon levelset function
        !   lambda          : scaling parameter

        ! Fields
        class(PolygonLevelsetFunction2DUDT), allocatable    :: targetplf
        real(R8)                                            :: lambda

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionPLF

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionPLF

    end type

    ! Face angle cost function
    type, extends(CostfunctionSOUDT) :: CostfunctionSOFAUDT

        ! Description
        !============
        ! This cost function penalizes the angle that vessel
        ! face make with respect to each other (actually the angle minus
        ! pi since that would give a straight line). This is done using
        ! an L2-norm. To this end, the vertex pairs of each edge are
        ! stored in an nvp-by-4 array. This allows to also penalize 
        ! non-subsequent edge pairs (though likely not necessary). The
        ! following fields are stored
        !
        !   nvpairs         : total number of vertex pairs
        !   vpairs          : vertex pairs (vessel vertex IDs)
        !   lambda          : scaling parameter
        !   wt              : weight vector to determine relative weight
        !                   of each vertex pair in the cost function

        ! Fields
        integer(I8), allocatable                :: vpairs(:, :)
        integer                                 :: nvpairs 
        real(R8)                                :: lambda 
        real(R8), allocatable                   :: wt(:)

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionSOFA

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionSOFA

        ! Data writing
        procedure :: WriteData              => WriteCostFunctionDataSOFA

    end type

    ! Cost function with all possible contributions (but vessel based only)
    type, extends(CostfunctionSOUDT) :: CostfunctionGeneralSOUDT

        ! Description
        !============
        ! Cost function that accounts for all possible combinations of 
        ! length ratio(s), angles, differences, ... The inclusion of a 
        ! cost function value is determined based on the value of the 
        ! scaling coefficient lambda. If this is zero or negative, the 
        ! contribution is not included. One should beware that if the 
        ! lambda values are not properly set in the input file, 
        ! contributions may be unexpectedly included since the default
        ! value for these contributions is non-zero. If no contributions
        ! would be included, the system is likely underdetermined, 
        ! leading to NaNs/divergence of the solver. 

        ! Fields
        type(CostfunctionPLFUDT)        :: cfv_plf
        type(CostFunctionSOFAUDT)       :: cfv_fa

        ! Switches
        logical                         :: doFA, doPLF

    contains 

        ! Initialization
        procedure :: Initialize         => InitializeCostfunctionGeneralSO

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionGeneralSO


    end type

    ! General goat-reduced cost function
    type, extends(CostfunctionSOUDT) :: CostFunctionGRUDT

        ! Description 
        !============
        ! General cost function for goat-reduced formulation. It has an
        ! allocatable cost function field which can take any of the 
        ! regular cost function options and the goat solver (the
        ! optimization problem is already saved as state variable 
        ! anyway). This allows for a nested formulation without 
        ! having to change solvers etc (i.e. KKT solver can be reused or
        ! any other if desired for the shape optimization problem). 
        ! Additionally, the full goat engine is added. Since the goat
        ! problem is already passed through to all the routines, the
        ! problem instance of the engine is simply equated to the 
        ! goat problem. 

        ! Note: one should *not* impose the goat constraints in 
        ! addition to the optimization problem 

        ! Note: since the hessian of the reduced problem is typically
        ! hard to compute, we use a hessian estimator for the reduced 
        ! cost function  (except for those contributions that can be 
        ! computed directly). This hessian estimator 

        ! Fields
        class(CostFunctionSOUDT), allocatable       :: costfunction 
        type(OptimizationEngineGDUDT)               :: goatengine 

    contains 

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionGR

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionGR

    end type

    
    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Abstract interfaces
    !====================
    ! Cost function
    abstract interface

        ! Cost function initialization
        subroutine InitializeCostfunctionSOINT(costfunction, goat, options)

            ! Description
            !============
            ! This routine should initialize all additional parameters
            ! that are needed to evaluate the cost function (e.g. the 
            ! vertex indices where the cost function is defined).
            
            ! Import
            import :: CostfunctionSOUDT, OptimizationProblemGDUDT, &
                CostFunctionOptionsSOUDT

            ! Declare 
            class(CostfunctionSOUDT)        :: costfunction 
            type(OptimizationProblemGDUDT)  :: goat
            type(CostFunctionOptionsSOUDT)  :: options

        end subroutine

        ! Cost function evaluation
        subroutine EvaluateCostFunctionSOINT(costfunction, J, gradJ, &
            hessJ, goat, dogradient, dohessian, designvariables, &
            varin, valuesin, dJdvarin, dgradJdvarin)

            ! Description
            !============
            ! Main routine to evaluate the cost function and its 
            ! derivative and hessian w.r.t. design variables. 

            ! Import
            import :: CostfunctionSOUDT, MySparseUDT, R8, &
                OptimizationProblemGDUDT, DesignVariablesSOUDT
            
            ! Declare
            class(CostfunctionSOUDT)        :: costfunction 
            real(R8)                        :: J 
            real(R8), allocatable           :: gradJ(:)
            type(MySparseUDT)               :: hessJ 
            type(OptimizationProblemGDUDT)  :: goat
            logical                         :: dogradient, dohessian
            class(DesignVariablesSOUDT)     :: designvariables

            character(*), intent(in), optional  :: varin 
            real(R8), intent(in), optional      :: valuesin(:)
            real(R8), allocatable, optional     :: dJdvarin(:) 
            type(MySparseUDT), optional         :: dgradJdvarin

        end subroutine

    end interface

    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                        DUMMY COST FUNCTION                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionDummy(costfunction, goat, options)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionDummyUDT)         :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Nothing to do here
        costfunction%doremesh = .false.
        
    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionDummy(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)


        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionDummyUDT)     :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

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

        ! Set outputs to zero
        J = 0
        gradJ = 0
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)

        ! Other derivatives
        !==================
        ! Initialize
        allocate(dJdvar(size(values)))
        dJdvar = 0
        dgradJdvar = SpZeros(size(gradJ), size(values)) ! jacobian, not gradient

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                         LEVELSET FUNCTION                        !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionPLF(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! goat. The cost function is defined as:
        ! 
        !   cfv = lambda sum_i 0.5*V(xv(i), yv(i))**2,
        ! 
        ! where V is the polygon levelset function of the desired
        ! vessel location. This plf is set up by reading in another
        ! polygon structure from a file specified in the options and
        ! constructing a vessel from that file. 
        
        ! It is very important to realize that this cost function does 
        ! not map specific vessel vertices to a specific location! 
        ! So results may vary (e.g. vessel vertices may collide, or if 
        ! different closed structures exist, some vessel points may lie
        ! on one or the other)

        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPLFUDT)           :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Auxiliary 
        integer(I8)                         :: filespec 
        type(VesselOptionsUDT)              :: vesseloptions
        type(VesselUDT)                     :: newvessel, triangulationvessel

        ! Data
        
        ! Initialize
        !===========
        ! Set remeshing
        costfunction%doremesh = .false.

        ! Set the scaling constant
        costfunction%lambda = options%plf%lambda

        ! Set the vessel options loading path (typically goat filepath)
        vesseloptions%inputfilepath = options%plf%vesselinputfilepath 

        ! Initialize the vessel options
        call vesseloptions%Set()

        ! Set the new vessel filepath
        vesseloptions%filepath = options%plf%newvesselfilepath

        ! Set file specifier
        filespec = 10
 
        ! Set up the vessel polygon
        !==========================
        ! Read in the vessel
        call ReadVessel(filespec, newvessel, triangulationvessel, vesseloptions)

        ! Extract data
        call ExtractVesselData(newvessel, vesseloptions)

        ! Set PLF
        !========
        costfunction%targetplf = newvessel%plfvessel
        
    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionPLF(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, the target plf is evaluated at the vessel coordinates 
        ! (which may or may not be design variables...)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionPLFUDT)       :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: xv, yv, val, valxx, &
            valxy, valyy, dvaldx, dvaldy, d2valdx2, d2valdxdy, d2valdy2

        ! Loop 
        integer(I8)                             :: k 

        ! Initialize
        !===========
        ! Associate
        associate(&
            ps          => goat%environment%vessel%polygonset,    &
            plf         => costfunction%targetplf)

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

        ! Get vessel coordinates
        call ps%GetVertices(xv, yv)

        ! Initialize
        nv = size(xv)
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)

        ! Value
        !======
        ! Simply evaluate the levelset function 
        allocate(val(nv))
        call plf%Evaluate(xv, yv, 0, 0, val)

        ! Compute cost function contribution
        J = 0.5*costfunction%lambda*sum(val**2)

        ! Gradient & Hessian
        !===================
        ! Compute
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! Same gradient contributions for both design variables, 
            ! since vessel coordinates come first in design variable 
            ! vector and since there are no contributions in terms of 
            ! goat coordinates

            ! Gradient
            !---------
            if (dogradient) then 
                
                ! Evaluate derivatives
                allocate(dvaldx(nv), dvaldy(nv))
                call plf%Evaluate(xv, yv, 1, 0, dvaldx)
                call plf%Evaluate(xv, yv, 0, 1, dvaldy)

                ! Evaluate contributions
                gradJ(1:nv)         = val*dvaldx 
                gradJ(nv+1:2*nv)    = val*dvaldy

            end if 

            ! Hessian
            !--------
            if (dohessian) then 

                ! Check 
                if (.not. dogradient) then 

                    ! Evaluate derivatives
                    allocate(dvaldx(nv), dvaldy(nv))
                    call plf%Evaluate(xv, yv, 1, 0, dvaldx)
                    call plf%Evaluate(xv, yv, 0, 1, dvaldy)

                end if 

                ! Evaluate second order derivatives
                allocate(d2valdx2(nv), d2valdxdy(nv), d2valdy2(nv))
                call plf%Evaluate(xv, yv, 2, 0, d2valdx2)
                call plf%Evaluate(xv, yv, 1, 1, d2valdxdy)
                call plf%Evaluate(xv, yv, 0, 2, d2valdy2)

                ! Evaluate contributions
                valxx = val*d2valdx2 + dvaldx**2
                valxy = val*d2valdxdy + dvaldx*dvaldy
                valyy = val*d2valdy2 + dvaldy**2

                ! Set hessian
                row = [(k, k = 1, nv), (k, k = 1, nv), &
                    (k, k = nv+1, 2*nv), (k, k = nv+1, 2*nv)]
                col = [(k, k = 1, nv), (k, k = nv+1, 2*nv), &
                    (k, k = 1, nv), (k, k = nv+1, 2*nv)]
                val = [valxx, valxy, valxy, valyy]
                hessJ = ConstructMySparse(row, col, val, designvariables%nphi, designvariables%nphi)

            end if 



        case default 

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionPLF: gradient ' // & 
                ' and hessian not implemented for design variables: ' // &
                designvariables%type)

        end select

        ! Scale
        !------
        gradJ = gradJ*costfunction%lambda 
        hessJ = hessJ*costfunction%lambda

        ! Other derivatives
        !==================
        ! Initialize
        allocate(dJdvar(size(values)))
        dJdvar = 0
        dgradJdvar = SpZeros(size(gradJ), size(values)) ! jacobian, not gradient

        ! Currently no derivatives w.r.t. any other variables
        select case (var)

        case ('goatvariables', 'no')

        case default 

            call gdErrorHandler('EvaluateCostFunctionPLF: unknown ' // & 
                'variable for derivative calculation: ' // var)

        end select

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

        ! Assocation termination
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                            FACE ANGLE                            !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionSOFA(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! goat. The cost function is defined as:
        ! 
        !   cfv = lambda sum_i (theta_i - pi)**2,
        ! 
        ! where theta_i is the angle that is made between two vertex
        ! pairs of the vessel. Note that this cost function promotes
        ! straight surfaces and may not yield satisfactory results
        ! if used on its own. 

        ! Note: we assume that all vessel polygons are closed
        
        ! Modules

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionSOFAUDT)          :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Auxiliary 
        integer(I8)                                 :: nvpairs
        integer(I8), allocatable, dimension(:, :)   :: vpairs  
        real(R8), allocatable, dimension(:)         :: xv, yv 

        ! Loop
        
        ! Initialize
        !===========
        ! Set remeshing
        costfunction%doremesh = .false.

        ! Associate
        associate(&
            opt       => options%fa,                        &   
            ps        => goat%environment%vessel%polygonset &
            )

        ! Set the scaling constant
        costfunction%lambda = opt%lambda

        ! Determine vertex pairs
        call goat%environment%vessel%GetVesselVertexPairs(vpairs, &
            opt%structureIDs, opt%vertIDs)
        nvpairs = size(vpairs, 1)

        ! Add
        !====
        allocate(costfunction%vpairs(nvpairs, 4))
        costfunction%vpairs(:, 1:2) = vpairs(:, 1:2)
        costfunction%vpairs(:, 3:4) = vpairs(:, 2:3)
        costfunction%nvpairs = nvpairs 

        ! Determine weights
        !==================
        ! Just equal to one for now
        allocate(costfunction%wt(nvpairs))
        costfunction%wt = 1 ! default

        ! Write data
        !===========
        call ps%GetVertices(xv, yv)
        call costfunction%WriteData(xv, yv)

        ! Housekeeping
        !=============
        end associate
        
    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionSOFA(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, the target plf is evaluated at the vessel coordinates 
        ! (which may or may not be design variables...)

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionSOFAUDT)      :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                             :: nv 
        integer(I8), allocatable, dimension(:)  :: row, col
        real(R8), allocatable, dimension(:)     :: x1v, x2v, x3v, x4v, &
            y1v, y2v, y3v, y4v, dx1v, dx2v, dy1v, dy2v, dpv, cpv, ratv, &
            thetav, xv, yv, valxx, valxy, valyx, valyy

        ! Loop 
        integer(I8)                             :: i, k 

        ! Initialize
        !===========
        ! Associate
        associate(&
            ps          => goat%environment%vessel%polygonset,    &
            vpairs      => costfunction%vpairs,     &
            nvpairs     => costfunction%nvpairs,    &
            lambda      => costfunction%lambda,     &
            wt          => costfunction%wt)

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

        ! Check design variables
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! All good

        case default

            ! All bad
            call gdErrorHandler('EvaluateCostFunctionSOFA: design variable ' // & 
                'type: "' // designvariables%type // '" not implemented')

        end select
                
        ! Get vessel coordinates
        call ps%GetVertices(xv, yv)

        ! Initialize
        nv = size(xv)
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)
        gradJ = 0

        ! Value
        !======
        ! Compute angle
        x1v = xv(vpairs(:, 1))
        x2v = xv(vpairs(:, 2))
        x3v = xv(vpairs(:, 3))
        x4v = xv(vpairs(:, 4))

        y1v = yv(vpairs(:, 1))
        y2v = yv(vpairs(:, 2))
        y3v = yv(vpairs(:, 3))
        y4v = yv(vpairs(:, 4))

        dx1v = x2v - x1v 
        dx2v = x4v - x3v 
        dy1v = y2v - y1v 
        dy2v = y4v - y3v 

        dpv = dx1v*dx2v + dy1v*dy2v 
        cpv = dx1v*dy2v - dx2v*dy1v 

        ratv = cpv/dpv 
        thetav = atan2(cpv, dpv)

        ! Compute cost function contribution
        J = 0.5*lambda*sum(thetav**2)

        ! Gradient & Hessian
        !===================
        ! Compute
        select case (designvariables%type)

        case ('vesselcoordinates', 'vesselcoordinates_goat')

            ! Same gradient contributions for both design variables, 
            ! since vessel coordinates come first in design variable 
            ! vector and since there are no contributions in terms of 
            ! goat coordinates

            ! Gradient
            !---------
            if (dogradient) then 

                ! Evaluate contributions
                do i = 1, nvpairs 
                    ! Unpack 
                    associate(&
                        dx1   => dx1v(i),       dx2     => dx2v(i),     &
                        dy1   => dy1v(i),       dy2     => dy2v(i),     &
                        dp    => dpv(i),        cp      => cpv(i),      &
                        rat   => ratv(i),       theta   => thetav(i),   &
                        wti   => wt(i)  &
                        )

                    gradJ(vpairs(i, 1)) = gradJ(vpairs(i, 1)) + &
                        -(theta*wti*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1) !x1
                    gradJ(vpairs(i, 2)) = gradJ(vpairs(i, 2)) + &
                        (theta*wti*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1) !x2
                    gradJ(vpairs(i, 3)) = gradJ(vpairs(i, 3)) + &
                        (theta*wti*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1) !x3
                    gradJ(vpairs(i, 4)) = gradJ(vpairs(i, 4)) + &
                        -(theta*wti*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1) !x4

                    gradJ(vpairs(i, 1)+nv) = gradJ(vpairs(i, 1)+nv) + &
                        (theta*wti*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1) !y1
                    gradJ(vpairs(i, 2)+nv) = gradJ(vpairs(i, 2)+nv) + &
                        -(theta*wti*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1) !y2
                    gradJ(vpairs(i, 3)+nv) = gradJ(vpairs(i, 3)+nv) + &
                        -(theta*wti*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1) !y3
                    gradJ(vpairs(i, 4)+nv) = gradJ(vpairs(i, 4)+nv) + &
                        (theta*wti*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1) !y4

                    end associate 
                end do 

            end if 

            ! Hessian
            !--------
            if (dohessian) then 

                ! Allocate
                allocate(valxx(nvpairs*16), valxy(nvpairs*16), &
                    valyx(nvpairs*16), valyy(nvpairs*16), &
                    row(nvpairs*16), col(nvpairs*16))

                ! Initialize
                k = 0

                ! Unpack 
                associate(&
                    dx1   => dx1v,       dx2     => dx2v,     &
                    dy1   => dy1v,       dy2     => dy2v,     &
                    dp    => dpv,        cp      => cpv,      &
                    rat   => ratv,       theta   => thetav,   &
                    wti   => wt  &
                    )

                ! v1v1
                row(k+1:k+nvpairs) = vpairs(:, 1)
                col(k+1:k+nvpairs) = vpairs(:, 1)
                valxx(k+1:k+nvpairs) =  (wti*(dy2/dp - (dx2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dx2**2*rat)/dp**2 - (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x1
                valxy(k+1:k+nvpairs) = (theta*wti*(dx2**2/dp**2 - dy2**2/dp**2 &
                    + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) - (wti*(dx2/dp &
                    + (dy2*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y1
                valyx(k+1:k+nvpairs) = valxy(k+1:k+nvpairs) ! y1x1
                valyy(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dy2**2*rat)/dp**2 + (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y1
                k = k + nvpairs

                ! v1v2
                row(k+1:k+nvpairs) = vpairs(:, 1)
                col(k+1:k+nvpairs) = vpairs(:, 2)
                valxx(k+1:k+nvpairs) =  - (wti*(dy2/dp &
                    - (dx2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dx2**2*rat)/dp**2 &
                    - (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x2
                valxy(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dx2**2/dp**2 &
                    - dy2**2/dp**2 + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y2
                valyx(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dx2**2/dp**2 &
                    - dy2**2/dp**2 + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1x2
                valyy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy2**2*rat)/dp**2 &
                    + (2*dx2*dy2)/dp**2))/(rat**2 + 1) - (wti*(dx2/dp &
                    + (dy2*rat)/dp)**2)/(rat**2 + 1)**2 !y1y2
                k = k + nvpairs

                ! v1v3
                row(k+1:k+nvpairs) = vpairs(:, 1)
                col(k+1:k+nvpairs) = vpairs(:, 3)
                valxx(k+1:k+nvpairs) = (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 &
                    - (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(rat/dp &
                    + (dx1*dy2)/dp**2 - (dx2*dy1)/dp**2 &
                    - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) !x1x3
                valxy(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y3
                valyx(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((dx1*dx2)/dp**2 &
                    - 1/dp + (dy1*dy2)/dp**2 + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1x3
                valyy(k+1:k+nvpairs) = - (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y3
                k = k + nvpairs

                ! v1v4
                row(k+1:k+nvpairs) = vpairs(:, 1)
                col(k+1:k+nvpairs) = vpairs(:, 4)
                valxx(k+1:k+nvpairs) =  (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1x4
                valxy(k+1:k+nvpairs) = - (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x1y4
                valyx(k+1:k+nvpairs) = (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((dx1*dx2)/dp**2 - 1/dp + (dy1*dy2)/dp**2 &
                    + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1x4
                valyy(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y1y4
                k = k + nvpairs

                ! v2v1
                row(k+1:k+nvpairs) = vpairs(:, 2)
                col(k+1:k+nvpairs) = vpairs(:, 1)
                valxx(k+1:k+nvpairs) = - (wti*(dy2/dp &
                    - (dx2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dx2**2*rat)/dp**2 &
                    - (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2x1
                valxy(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dx2**2/dp**2 &
                    - dy2**2/dp**2 + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2y1
                valyx(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dx2**2/dp**2 &
                    - dy2**2/dp**2 + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2x1
                valyy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy2**2*rat)/dp**2 + (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    - (wti*(dx2/dp + (dy2*rat)/dp)**2)/(rat**2 + 1)**2 !y2y1
                k = k + nvpairs

                ! v2v2
                row(k+1:k+nvpairs) = vpairs(:, 2)
                col(k+1:k+nvpairs) = vpairs(:, 2)
                valxx(k+1:k+nvpairs) =  (wti*(dy2/dp &
                    - (dx2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dx2**2*rat)/dp**2 - (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2x2
                valxy(k+1:k+nvpairs) = (theta*wti*(dx2**2/dp**2 - dy2**2/dp**2 &
                    + (2*dx2*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (wti*(dx2/dp + (dy2*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2y2
                valyx(k+1:k+nvpairs) = valxy(k+1:k+nvpairs)
                valyy(k+1:k+nvpairs) = (wti*(dx2/dp + (dy2*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dy2**2*rat)/dp**2 + (2*dx2*dy2)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2y2
                k = k + nvpairs

                ! v2v3
                row(k+1:k+nvpairs) = vpairs(:, 2)
                col(k+1:k+nvpairs) = vpairs(:, 3)
                valxx(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2x3
                valxy(k+1:k+nvpairs) = - (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2y3
                valyx(k+1:k+nvpairs) = (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((dx1*dx2)/dp**2 - 1/dp + (dy1*dy2)/dp**2 &
                    + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2x3
                valyy(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2y3
                k = k + nvpairs

                ! v2v4
                row(k+1:k+nvpairs) = vpairs(:, 2)
                col(k+1:k+nvpairs) = vpairs(:, 4)
                valxx(k+1:k+nvpairs) =  (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy2/dp - (dx2*rat)/dp))/(rat**2 + 1)**2 &
                    - (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(rat/dp &
                    + (dx1*dy2)/dp**2 - (dx2*dy1)/dp**2 &
                    - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) !x2x4
                valxy(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 !x2y4
                valyx(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((dx1*dx2)/dp**2 &
                    - 1/dp + (dy1*dy2)/dp**2 + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2x4
                valyy(k+1:k+nvpairs) = - (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !y2y4
                k = k + nvpairs

                ! v3v1
                row(k+1:k+nvpairs) = vpairs(:, 3)
                col(k+1:k+nvpairs) = vpairs(:, 1)
                valxx(k+1:k+nvpairs) =  - (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    - (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3x1
                valxy(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((dx1*dx2)/dp**2 &
                    - 1/dp + (dy1*dy2)/dp**2 + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3y1
                valyx(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3x1
                valyy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 &
                    - (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(rat/dp &
                    + (dx1*dy2)/dp**2 - (dx2*dy1)/dp**2 &
                    - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) !y3y1
                k = k + nvpairs

                ! v3v2
                row(k+1:k+nvpairs) = vpairs(:, 3)
                col(k+1:k+nvpairs) = vpairs(:, 2)
                valxx(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3x2
                valxy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((dx1*dx2)/dp**2 - 1/dp + (dy1*dy2)/dp**2 &
                    + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !x3y2
                valyx(k+1:k+nvpairs) = - (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3x2
                valyy(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3y2
                k = k + nvpairs

                ! v3v3
                row(k+1:k+nvpairs) = vpairs(:, 3)
                col(k+1:k+nvpairs) = vpairs(:, 3)
                valxx(k+1:k+nvpairs) =  (wti*(dy1/dp + (dx1*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dx1**2*rat)/dp**2 + (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3x3
                valxy(k+1:k+nvpairs) = (theta*wti*(dy1**2/dp**2 - dx1**2/dp**2 &
                    + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3y3
                valyx(k+1:k+nvpairs) = valxy(k+1:k+nvpairs)
                valyy(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dy1**2*rat)/dp**2 - (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3y3
                k = k + nvpairs

                ! v3v4
                row(k+1:k+nvpairs) = vpairs(:, 3)
                col(k+1:k+nvpairs) = vpairs(:, 4)
                valxx(k+1:k+nvpairs) =  (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dx1**2*rat)/dp**2 + (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    - (wti*(dy1/dp + (dx1*rat)/dp)**2)/(rat**2 + 1)**2 !x3x4
                valxy(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dy1**2/dp**2 &
                    - dx1**2/dp**2 + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x3y4
                valyx(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dy1**2/dp**2 &
                    - dx1**2/dp**2 + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3x4
                valyy(k+1:k+nvpairs) = - (wti*(dx1/dp - (dy1*rat)/dp)**2)/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy1**2*rat)/dp**2 - (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y3y4
                k = k + nvpairs

                ! v4v1
                row(k+1:k+nvpairs) = vpairs(:, 4)
                col(k+1:k+nvpairs) = vpairs(:, 1)
                valxx(k+1:k+nvpairs) =  (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4x1
                valxy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((dx1*dx2)/dp**2 - 1/dp + (dy1*dy2)/dp**2 &
                    + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx2/dp + (dy2*rat)/dp))/(rat**2 + 1)**2 !x4y1
                valyx(k+1:k+nvpairs) = - (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4x1
                valyy(k+1:k+nvpairs) = (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    + (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4y1
                k = k + nvpairs

                ! v4v2
                row(k+1:k+nvpairs) = vpairs(:, 4)
                col(k+1:k+nvpairs) = vpairs(:, 2)
                valxx(k+1:k+nvpairs) =  - (theta*wti*(rat/dp + (dx1*dy2)/dp**2 &
                    - (dx2*dy1)/dp**2 - (2*dx1*dx2*rat)/dp**2))/(rat**2 + 1) &
                    - (wti*(dy1/dp + (dx1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*((2*dx2*rat**2)/dp &
                    - (2*dy2*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4x2
                valxy(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*((dx1*dx2)/dp**2 &
                    - 1/dp + (dy1*dy2)/dp**2 + (2*dx1*dy2*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dy2*rat**2)/dp + (2*dx2*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4y2
                valyx(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)*(dy2/dp &
                    - (dx2*rat)/dp))/(rat**2 + 1)**2 + (theta*wti*(1/dp &
                    - (dx1*dx2)/dp**2 - (dy1*dy2)/dp**2 &
                    + (2*dx2*dy1*rat)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dx2*rat**2)/dp - (2*dy2*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4x2
                valyy(k+1:k+nvpairs) = (theta*wti*((2*dy2*rat**2)/dp &
                    + (2*dx2*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 &
                    - (wti*(dx1/dp - (dy1*rat)/dp)*(dx2/dp &
                    + (dy2*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(rat/dp &
                    + (dx1*dy2)/dp**2 - (dx2*dy1)/dp**2 &
                    - (2*dy1*dy2*rat)/dp**2))/(rat**2 + 1) !y4y2
                k = k + nvpairs

                ! v4v3
                row(k+1:k+nvpairs) = vpairs(:, 4)
                col(k+1:k+nvpairs) = vpairs(:, 3)
                valxx(k+1:k+nvpairs) = (theta*wti*((2*dx1*rat**2)/dp &
                    + (2*dy1*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dx1**2*rat)/dp**2 + (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    - (wti*(dy1/dp + (dx1*rat)/dp)**2)/(rat**2 + 1)**2 !x4x3
                valxy(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dy1**2/dp**2 - dx1**2/dp**2 &
                    + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) + (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dy1/dp + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4y3
                valyx(k+1:k+nvpairs) = (wti*(dy1/dp + (dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 - (theta*wti*(dy1**2/dp**2 &
                    - dx1**2/dp**2 + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4x3
                valyy(k+1:k+nvpairs) = - (wti*(dx1/dp &
                    - (dy1*rat)/dp)**2)/(rat**2 + 1)**2 - (theta*wti*((2*dy1**2*rat)/dp**2 &
                    - (2*dx1*dy1)/dp**2))/(rat**2 + 1) - (theta*wti*((2*dy1*rat**2)/dp &
                    - (2*dx1*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4y3
                k = k + nvpairs

                ! v4v4
                row(k+1:k+nvpairs) = vpairs(:, 4)
                col(k+1:k+nvpairs) = vpairs(:, 4)
                valxx(k+1:k+nvpairs) =  (wti*(dy1/dp + (dx1*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dx1**2*rat)/dp**2 + (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    - (theta*wti*((2*dx1*rat**2)/dp + (2*dy1*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4x4
                valxy(k+1:k+nvpairs) = (theta*wti*(dy1**2/dp**2 - dx1**2/dp**2 &
                    + (2*dx1*dy1*rat)/dp**2))/(rat**2 + 1) - (wti*(dy1/dp &
                    + (dx1*rat)/dp)*(dx1/dp - (dy1*rat)/dp))/(rat**2 + 1)**2 &
                    - (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dy1/dp &
                    + (dx1*rat)/dp))/(rat**2 + 1)**2 !x4y4
                valyx(k+1:k+nvpairs) = valxy(k+1:k+nvpairs)
                valyy(k+1:k+nvpairs) = (wti*(dx1/dp - (dy1*rat)/dp)**2)/(rat**2 + 1)**2 &
                    + (theta*wti*((2*dy1**2*rat)/dp**2 - (2*dx1*dy1)/dp**2))/(rat**2 + 1) &
                    + (theta*wti*((2*dy1*rat**2)/dp - (2*dx1*rat)/dp)*(dx1/dp &
                    - (dy1*rat)/dp))/(rat**2 + 1)**2 !y4y4
                k = k + nvpairs

                ! Construct
                hessJ = ConstructMySparse([row, row, row + nv, row + nv], &
                    [col, col + nv, col, col + nv], [valxx, valxy, valyx, valyy], &
                    designvariables%nphi, designvariables%nphi)

                ! Housekeeping
                end associate

            end if 



        case default 

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionPLF: gradient ' // & 
                ' and hessian not implemented for design variables: ' // &
                designvariables%type)

        end select

        ! Scale
        !------
        gradJ = gradJ*costfunction%lambda 
        hessJ = hessJ*costfunction%lambda

        ! Other derivatives
        !==================
        ! Initialize
        allocate(dJdvar(size(values)))
        dJdvar = 0
        dgradJdvar = SpZeros(size(gradJ), size(values)) ! jacobian, not gradient

        ! Currently no derivatives w.r.t. any other variables
        select case (var)

        case ('goatvariables', 'no')

        case default 

            call gdErrorHandler('EvaluateCostFunctionPLF: unknown ' // & 
                'variable for derivative calculation: ' // var)

        end select

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

        ! Assocation termination
        end associate

    end subroutine

    ! Cost function data writing
    subroutine WriteCostFunctionDataSOFA(costfunction, xv, yv)

        ! Description
        !============
        ! This routine writes out the cost function data (vertex pairs)
        ! for post-processing

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionSOFAUDT)          :: costfunction 
        real(R8), intent(in)                :: xv(:), yv(:)

        ! Auxiliary
        integer(I8)                     :: ncol, nrow 

        integer(I8), allocatable        :: IDn(:, :) 
        real(R8), allocatable           :: xn(:, :), yn(:, :)
        character(:), allocatable       :: filename 

        ! Loop
        integer(I8)                     :: j 

        ! Initialize
        !===========
        ! Set filename
        filename = 'so_cfv_sofa_vpairs'

        ! Allocate
        nrow = size(costfunction%vpairs, 1)
        ncol = size(costfunction%vpairs, 2)
        allocate(IDn(nrow, ncol), xn(nrow, ncol), yn(nrow, ncol))

        ! Unpack
        associate(vpairs      => costfunction%vpairs)

        ! Loop
        xn = 0
        yn = 0
        do j = 1, nrow 
            IDn(j, :) = vpairs(j, :) 
            xn(j, :) = xv(vpairs(j, :))
            yn(j, :) = yv(vpairs(j, :))
        end do

        ! Call writer
        !============
        call WriteVertexPairData(IDn, xn, yn, filename)

        ! Housekeeping
        !=============
        end associate
        deallocate(IDn, xn, yn)

    end subroutine

    !------------------------------------------------------------------!
    !                             GENERAL                              !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGeneralSO(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the cost function and its parameters based on the 
        ! grid, magnetic field, and environment structures. 

        ! Simply call the initialization of the original lenght ratio
        ! cost function. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralSOUDT)     :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options
        
        ! Initialize
        !===========
        ! Set remesh
        costfunction%doremesh = .false.

        ! Set evaluation switches
        costfunction%doPLF      = .false.
        costfunction%doFA       = .false.

        ! Check based on cost function type
        select case (costfunction%type)

        case ('general')

            ! Include all contributions

        case default 

            ! Throw error
            call gdErrorHandler('Unknown cost function type')

        end select

        ! Initialize if necessary
        if (options%PLF%lambda > 0) then
            costfunction%doPLF = .true.
            call costfunction%cfv_plf%Initialize(goat, options)
        end if 
        if (options%FA%lambda > 0) then 
            costfunction%doFA = .true.
            call costfunction%cfv_fa%Initialize(goat, options)
        end if 

    end subroutine

    ! Cost function evaluation
    subroutine EvaluateCostFunctionGeneralSO(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, &
        designvariables, varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! Here, we simply call the same cost function twice, but switch
        ! the order of the indices and recompute the bias. 

        ! Notes:
        !=======
        ! Possible future performance improvements:
        ! - Allocating hessian stuff only once and storing indices, 
        ! since they don't change
        ! - Instead of recomputing auxiliary variables, store them. May
        ! not actually be better in terms of computational time, but 
        ! may lead to shorter and hence better maintainable code. 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGeneralSOUDT) :: costfunction
        real(R8)                        :: J, Jtemp
        real(R8), allocatable           :: gradJ(:), gradJtemp(:)
        type(MySparseUDT)               :: hessJ, hessJtemp
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:), dJdvartemp(:) 
        type(MySparseUDT)                   :: dgradJdvar, dgradJdvartemp

        ! Loop variables

        ! Auxiliary
                                        
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

        ! Cost function
        J = 0
        Jtemp = 0

        ! Gradient
        gradJ = 0
        allocate(gradJtemp(size(gradJ)), dJdvartemp(size(values)), &
            dJdvar(size(values)))
        gradJtemp = 0
        dJdvartemp = 0
        dJdvar = 0

        ! Hessian
        hessJtemp%nrow = hessJ%nrow 
        hessJtemp%ncol = hessJ%ncol 
        dgradJdvar = SpZeros(size(gradJ), size(values))
        dgradJdvartemp = dgradJdvar
        
        ! Allocate initially to avoid errors 
        hessJ = SpZeros(designvariables%nphi, designvariables%nphi)
            
        ! Compute cost function
        !======================
        ! PLF
        if (costfunction%doPLF) then 
            ! Compute
            call costfunction%cfv_plf%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, goat, dogradient, &
                dohessian, designvariables, var, values, dJdvartemp, dgradJdvartemp)
            
            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp
            dJdvar  = dJdvar + dJdvartemp 
            dgradJdvar  = dgradJdvar + dgradJdvartemp

            if (any(.not. ieee_is_finite(gradJtemp))) then 
                print *, 'Non-finite values in gradJ for PLF'
            end if

            ! Deallocate
            call hessJtemp%Deallocate()
            call dgradJdvartemp%Deallocate()
        end if 
        
        ! Face angle 
        if (costfunction%doFA) then 
            ! Compute
            call costfunction%cfv_fa%Evaluate(Jtemp, gradJtemp, &
                hessJtemp, goat, dogradient, &
                dohessian, designvariables, var, values, dJdvartemp, dgradJdvartemp)

            ! Add
            J       = J + Jtemp 
            gradJ   = gradJ + gradJtemp
            hessJ   = hessJ + hessJtemp
            dJdvar  = dJdvar + dJdvartemp 
            dgradJdvar  = dgradJdvar + dgradJdvartemp

            if (any(.not. ieee_is_finite(gradJtemp))) then 
                print *, 'Non-finite values in gradJ for FA'
            end if

            ! Deallocate
            call hessJtemp%Deallocate()
            call dgradJdvartemp%Deallocate()
        end if 

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if

    end subroutine

    !------------------------------------------------------------------!
    !                     GOAT REDUCED COST FUNCTION                   !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGR(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the goat reduced cost function. This basically means
        ! initializing the goat engine and the costfunction type. The 
        ! cost function is then allocated to the correct type and 
        ! initialized here as well. 

        ! Arguments
        class(CostfunctionGRUDT)            :: costfunction
        type(OptimizationProblemGDUDT)      :: goat
        type(CostFunctionOptionsSOUDT)      :: options

        ! Initialize goat engine
        !=======================
        ! Set the input filepaths for the driver
        costfunction%goatengine%inputfilepath = goat%inputfilepath
        costfunction%goatengine%inputfileprefix = 'gd.'
        
        ! Initialize the driver
        call costfunction%goatengine%SetupOptimizationDriver()

        ! Initialize the problem
        costfunction%goatengine%problem = goat

        ! Initialize the solver
        call costfunction%goatengine%solver%Initialize()

        ! Initialize costfunction
        !========================
        ! Set remesh
        costfunction%doremesh = .false.

        ! Check which cost function it is
        select case (options%type)

            case ('PLF')

                ! Allocate 
                allocate(CostFunctionPLFUDT::costfunction%costfunction)

                ! Set type
                costfunction%costfunction%type = 'PLF'

            case ('SOFA')

                ! Allocate
                allocate(CostFunctionSOFAUDT::costfunction%costfunction)

                ! Set type
                costfunction%costfunction%type = 'SOFA'

            case ('no')

                ! Zero cost function
                allocate(CostFunctionDummyUDT::costfunction%costfunction)

                ! Set type
                costfunction%costfunction%type = 'no'

            case ('general')
                
                ! General cost function
                allocate(CostFunctionGeneralSOUDT::costfunction%costfunction)

                ! Set type
                costfunction%costfunction%type = 'general'

        case default 

            ! Throw error
            call gdErrorHandler('InitializeCostFunctionGR: unknown cost' // & 
                'function type: ' // options%type)

        end select

        ! Initialize the cost function further
        call costfunction%costfunction%Initialize(goat, options)

    end subroutine

    ! Evaluation
    subroutine EvaluateCostFunctionGR(costfunction, J, gradJ, hessJ, &
        goat, dogradient, dohessian, designvariables, &
        varin, valuesin, dJdvarin, dgradJdvarin)

        ! Description
        !============
        ! Evaluate the cost function, the gradient and its hessian. 
        ! First, the goat equations are solved by calling the driver (
        ! it is assumed that goat is up to date). Afterwards, the cost
        ! function is evaluated and its gradient is computed by 
        ! applying a discrete adjoint approach to account for the
        ! goat constraints. The contributions of the goat constraints to
        ! the hessian of the problem are not accounted for... 

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGRUDT)        :: costfunction 
        real(R8)                        :: J
        real(R8), allocatable           :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian 
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin 
        real(R8), intent(in), optional      :: valuesin(:)
        real(R8), allocatable, optional     :: dJdvarin(:) 
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8), allocatable               :: values(:)
        real(R8), allocatable               :: dJdvar(:) 
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8)                                 :: flag, errstat
        real(R8)                                    :: goatcfv 

        real(R8), allocatable, dimension(:)         :: goatvariables, &
            gradJR, gradJgoat, lambdaG 

        type(MySparseUDT)                           :: hessJR, &
            jacGgoat, jacGdes, gradGdes, gradGgoat

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

        ! Initialize others
        allocate(gradJR(designvariables%nphi))
        gradJR = 0

        ! Construct goat variables
        goatvariables = [goat%designvariables%phi, goat%lambda, goat%mu]

        ! Initialize others
        allocate(gradJgoat(size(goatvariables)))

        ! Solve goat 
        !===========
        associate(goatproblem       => costfunction%goatengine%problem)
        select type (goatproblem)

        type is (OptimizationProblemGDUDT)

            ! Update the problem with the latest goat
            goatproblem = goat 

            ! Call the driver
            call ErrorStack%StartTrack()
            call costfunction%goatengine%solver%SolveOptimizationProblem(goatproblem)

            ! Check for any issues, if so - don't update the problem!
            errstat = ErrorStack%ErrorState()
            call ErrorStack%EndTrack()
            if (errstat > 0) then 
                ! No need to recall error, but need to check how to adjust output
                goatcfv = goatproblem%monitor%convnorm(goatproblem%monitor%itopt)
                call costfunction%goatengine%solver%Initialize()

                ! Raise the remesh flag
                costfunction%doremesh = .true. 
            else
                ! Update goat
                goat = goatproblem 
                goatcfv = goatproblem%monitor%convnorm(goatproblem%monitor%itopt)
            end if 

        class default 

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionGR: unexpected goat type')

        end select
        end associate

        ! Compute reduced cost function
        !==============================
        ! Call evaluation routine
        hessJR = SpZeros(designvariables%nphi, designvariables%nphi)
        call costfunction%costfunction%Evaluate(J, gradJR, hessJR, &
            goat, dogradient, dohessian, designvariables, &
            'goatvariables', goatvariables, gradJgoat)

        ! Compute gradient 
        !=================
        ! Evaluate goat jacobian w.r.t. goat variables
        call goat%EvaluateJacobian('goatvariables', goatvariables, jacGgoat)

        ! Evaluate goat jacobian w.r.t. design variables
        select case (designvariables%type)

        case ('vesselcoordinates')

            ! Derivatives w.r.t. vessel coordinate
            call goat%EvaluateJacobian('vesselcoordinates', &
                designvariables%phi, jacGdes)

        case ('vesselcoordinates_goat')

            ! Illegal, throw error
            call gdErrorHandler('EvaluteCostFunctionGR: goat variables' // & 
                ' cannot be present as explicit design variables, ' // & 
                ' not supported')

        case default 

            ! Unknown 
            call gdErrorHandler('EvaluateCostFunctionGR: design variable' // & 
                ' gradient w.r.t. design variables of type: ' // & 
                designvariables%type // ' are not implemented')

        end select

        ! Add goat residual to cost function
        J = J + goatcfv

        ! Compute lagrange multipliers
        gradGgoat = jacGgoat%Transpose()
        call SolveSparseLinearSystemDI(gradGgoat, -gradJgoat, lambdaG, flag)

        ! Compute gradient
        gradGdes = jacGdes%Transpose()
        gradJ = gradJR + gradGdes%MatrixVectorProduct(lambdaG)

        ! Compute hessian
        !================
        ! Only take the raw cost function contribution as hessian - 
        ! perhaps replace by hessian estimator in the future? 
        hessJ = hessJR

        ! Housekeeping
        !=============
        ! Optional arguments
        if (present(dJdvarin)) then 
            dJdvarin = dJdvar 
        end if 
        if (present(dgradJdvarin)) then 
            dgradJdvarin = dgradJdvar
        end if




    end subroutine

end module
