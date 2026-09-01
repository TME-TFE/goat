!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains the cost function implementation for shape
! optimization with SOLPS. It assumes that (at least) the b25 library
! with cost function (gradient) evaluation routines is available and
! linked to. It also makes use of the types defined in b2mod.


module sosmod_costfunction

    ! Initialize
    !============
    ! Load modules
    use b2mod_switches_diff
    use b2us_map_diff
    use b2mod_ad_diff, only : nncf
    use b2mod_main_diff, only : b2mn_init_b
    use b2us_geo_diff
    use b2us_plasma_diff
    !use b2mod_agdr_diff, only : b2agdr_opt, b2agdr_opt_b, b2agdr_init, &
    !    b2agdr_init_b, b2agdr_write, b2agdr_fin, b2agdr_fin_b
    use b2mod_par_opt_diff
    use somod_costfunction
    use gdmod_optimizationengine
    use b2mod_costfunction
    use b2mod_state_constr_diff
    use b2mod_costfunction_data

    ! Rename precision ...
    use b2mod_types, only: R8_B25 => R8
    use mod_precision, only: R8_G => R8, I8_G => I8





    ! The usual
    implicit none
    save
    public


    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    ! Derived types
    !==============
    ! Basic solps cost function (not inherited, so has to be included
    ! in another cost function)
    type :: CostFunctionSOLPSUDT

        ! Description
        !============
        ! This cost function is a reduced cost function that is
        ! fully based on a cost function definition from SOLPS. It is
        ! assumed that there is a routine called
        ! 'EvaluateCostFunctionGradient' that returns the cost
        ! function value and its gradient. Our routines wrap around it
        ! and this UDT defines the required data structures from SOLPS
        ! to initialize, evaluate and update the cost function.

        ! Note: we follow the same naming convention as defned in the
        ! SOLPS routines for (hopefully) clarity

        ! Fields
        type(geometry)          :: g ! geometry and grid
        type(geometry_diff)     :: gb ! geometry derivatives
        type(b2state)           :: st ! state
        type(b2state_diff)      :: stb ! state derivatives
        type(b2stateext)        :: state_ext
        type(b2stateext_diff)   :: state_extb
        type(b2average)         :: state_avg
        type(b2average)         :: state_avgb
        type(switches)          :: switch
        type(switches_diff)     :: switchb
        type(mapping)           :: mpg
        type(mapping_diff)      :: mpgb
        type(state_constraint)  :: stc
        type(state_constraint_diff)     :: stcb
        type(costfunction_data)     :: cf_data

    contains

        ! Initialization
        procedure :: Initialize             => InitializeCostfunctionSOLPS

        ! Evaluation
        procedure :: Evaluate               => EvaluateCostFunctionSOLPS

    end type

    ! Goat and SOLPS reduced cost function
    type, extends(CostFunctionGRUDT)    :: CostFunctionGSRUDT

        ! Description
        !============
        ! This cost function type inherits from the goat reduced
        ! optimization cost function and should be used to incorporate
        ! the SOLPS cost function together with the goat. Additional
        ! cost function contributions from the purely goat side may be
        ! added in the future. Here, we need to overwrite the GR
        ! initialize and evaluation routines.

        ! Additional fields
        type(CostFunctionSOLPSUDT)      :: cfvsolps

    contains

        ! Initialization
        procedure :: Initialize         => InitializeCostFunctionGSR

        ! Evaluation
        procedure :: Evaluate           => EvaluateCostFunctionGSR


    end type


    contains

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                        SOLPS COST FUNCTION                       !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionSOLPS(costfunction)

        ! Description
        !============
        ! Initialize the SOLPS cost function. This is simply a wrapper
        ! for the InitializeCostFunction.F90 routine declared on the
        ! SOLPS side.

        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionSOLPSUDT)         :: costfunction

        ! Call initializer
        !=================
        call InitializeCostFunction(costfunction%switch, &
            costfunction%switchb, costfunction%g, costfunction%gb, &
            costfunction%mpg, costfunction%mpgb, costfunction%st, &
            costfunction%stb, costfunction%state_ext, costfunction%state_extb, &
            costfunction%state_avg, costfunction%state_avgb, &
            costfunction%stc, costfunction%stcb, costfunction%cf_data)

    end subroutine

     ! Evaluation
    subroutine EvaluateCostFunctionSOLPS(costfunction, J, gradJ, dogradient)

        ! Description
        !============
        ! Evaluate the SOLPS cost function and its gradient (hessian is
        ! not returned). The gradient is assumed to be with respect to
        ! the grid coordinates (vxx, vxy), the psi values at the
        ! vertices (vxfpsi), the derivatives of the psi values at the
        ! vertices (vxbx, vxby) and the ffbz value. This routine is a
        ! wrapper for the EvaluateCostFunctionGradient routine

        ! Use modules
        !============


        ! Declare variables
        !==================
        ! Arguments
        class(CostFunctionSOLPSUDT)                 :: costfunction
        real(R8_G), intent(out)                  :: J
        real(R8_G), allocatable, intent(out)     :: gradJ(:)
        logical, intent(in)                      :: dogradient

        ! Auxiliary
        real(R8_B25)                        :: J1(nncf), J1b(nncf)


        ! Initialize
        !===========
        ! Associate
        associate(&
            switch          => costfunction%switch,         &
            switchb         => costfunction%switchb,        &
            g               => costfunction%g,              &
            gb              => costfunction%gb,             &
            mpg             => costfunction%mpg,            &
            mpgb            => costfunction%mpgb,           &
            st              => costfunction%st,             &
            state_ext       => costfunction%state_ext,      &
            state_extb      => costfunction%state_extb,     &
            state_avg       => costfunction%state_avg,      &
            state_avgb      => costfunction%state_avgb,     &
            stb             => costfunction%stb,            &
            stc             => costfunction%stc,            &
            stcb            => costfunction%stcb,           &
            cf_data         => costfunction%cf_data         &
            )

        ! Evaluate
        !=========
        ! Call routine
        if (dogradient) then
            ! Evaluate the gradient as well
            call EvaluateCostFunctionGradient(switch, switchb, g, gb, &
                mpg, mpgb, st, stb, state_ext, state_extb, state_avg, &
                state_avgb, J1, J1b, stc, stcb, cf_data)

            ! Extract and cast into our precision format
            J = real(J1(1), kind=R8_G) ! assumed first entry is total cost function
            gradJ = real([gb%vxx, gb%vxy, gb%vxfpsi, gb%vxbx, gb%vxby, &
                gb%vxffbz], kind=R8_G)
        else
            ! Evaluate only the cost function
            call EvaluateCostfunction(switch, g, mpg, st, &
                state_ext, state_avg, J1, stc, cf_data)

            ! Extract and cast into our precision format
            J = real(J1(1), kind=R8_G) ! assumed first entry is total cost function
            gradJ = 0*real([g%vxx, g%vxy, g%vxfpsi, g%vxbx, g%vxby, &
                g%vxffbz], kind=R8_G)
        end if

        ! Housekeeping
        !=============
        end associate

    end subroutine

    !------------------------------------------------------------------!
    !                GOAT AND SOLPS REDUCED COST FUNCTION              !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializeCostFunctionGSR(costfunction, goat, options)

        ! Description
        !============
        ! Initialize the goat reduced cost function. This basically means
        ! initializing the goat engine and the costfunction type. The
        ! cost function is then allocated to the correct type and
        ! initialized here as well.

        ! Arguments
        class(CostfunctionGSRUDT)           :: costfunction
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

        ! Initialize GOAT-based costfunction
        !===================================
        ! Check which cost function it is
        select case (options%type)

        case ('PLF')

            ! Allocate
            allocate(CostfunctionPLFUDT::costfunction%costfunction)

            ! Set type
            costfunction%costfunction%type = 'PLF'

        case ('SOFA')

            ! Allocate
            allocate(CostFunctionSOFAUDT::costfunction%costfunction)

            ! Set type
            costfunction%costfunction%type = 'SOFA'

        case ('general')

            ! General cost function
            allocate(CostFunctionGeneralSOUDT::costfunction%costfunction)

            ! Set type
            costfunction%costfunction%type = 'general'

        case default

            ! Throw error
            call gdErrorHandler('InitializeCostFunctionGRS: unknown cost' // &
                'function type: ' // options%type)

        end select

        ! Initialize the cost function further
        call costfunction%costfunction%Initialize(goat, options)

        ! Initialize SOLPS cost function
        !===============================
        ! Simply call initializer
        call costfunction%cfvsolps%Initialize()

    end subroutine

    ! Evaluation
    subroutine EvaluateCostFunctionGSR(costfunction, J, gradJ, hessJ, &
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

        ! Note: for the gradient, we need to account for the dependency
        ! of the psi values (and first order derivatives) on the
        ! grid coordinates in solps. The partial derivatives are
        ! available through the differentiated structures (geo_diff
        ! types for example), but we need to differentiate again w.r.t.
        ! the coordinates here. Additionally, we have to update the
        ! coordinate-dependent structures before evaluating the solps
        ! cost function. Note that we also update the grid cell
        ! coordinates, even though they are recomputed afterwards in
        ! solps - reason is that there is a vertex ordening step that
        ! (for solps reasons) comes before recomputation of the cell
        ! centers, which requires an estimate of the cell center
        ! coordinates. This also means that there is no gradient
        ! contribution of the cell center coordinates that has to be
        ! accounted for.

        ! Declare variables
        !==================
        ! Arguments
        class(CostfunctionGSRUDT)       :: costfunction
        real(R8_G)                      :: J
        real(R8_G), allocatable         :: gradJ(:) ! assumed initialized
        type(MySparseUDT)               :: hessJ ! assumed in initialized
        type(OptimizationProblemGDUDT)  :: goat
        logical                         :: dogradient, dohessian
        class(DesignVariablesSOUDT)     :: designvariables

        ! Optional arguments
        character(*), intent(in), optional  :: varin
        real(R8_G), intent(in), optional    :: valuesin(:)
        real(R8_G), allocatable, optional   :: dJdvarin(:)
        type(MySparseUDT), optional         :: dgradJdvarin

        character(:), allocatable           :: var
        real(R8_G), allocatable             :: values(:)
        real(R8_G), allocatable             :: dJdvar(:)
        type(MySparseUDT)                   :: dgradJdvar

        ! Auxiliary
        integer(I8_G)                       :: flag
        integer(I8_G), allocatable          :: tv(:)

        real(R8_G)                                  :: Js, Jg
        real(R8_G), allocatable, dimension(:)       :: goatvariables, &
            gradJg, gradJgoat, lambdaG, gradJs, dpsidx, dpsidy, &
            d2psidx2, d2psidxdy, d2psidy2, gradJsolps, psi, gradR

        type(MySparseUDT)                           :: hessJg, &
            jacGgoat, jacGdes, gradGdes, gradGgoat

        ! Loop
        integer(I8_G)                       :: i

        ! Error handling
        integer                             :: errstat

        ! Initialize
        !===========
        ! Associate for ease
        associate(nv            => goat%grid%vert%ntot, &
            MFinterp            => goat%magneticField%interp)

        ! Allocate
        allocate(psi(nv), dpsidx(nv), dpsidy(nv), d2psidx2(nv), d2psidy2(nv), &
            d2psidxdy(nv))

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
        allocate(gradJg(designvariables%nphi))
        gradJg = 0

        ! Construct goat variables
        goatvariables = [goat%designvariables%phi, goat%lambda, goat%mu]

        ! Initialize others
        allocate(gradJgoat(size(goatvariables)))

        ! Solve goat
        !===========


        ! Try to solve
        associate(goatproblem       => costfunction%goatengine%problem)
        select type (goatproblem)

        type is (OptimizationProblemGDUDT)

            ! Update the problem with the latest goat
            goatproblem = goat

            ! Keep track of errors
            call ErrorStack%StartTrack()

            ! Call the driver
            call costfunction%goatengine%solver%SolveOptimizationProblem(goatproblem)

            ! Check if an error was found, if so: call error (softly) and
            ! set cost function value and gradient to inf
            errstat = ErrorStack%ErrorState()
            call ErrorStack%EndTrack()
            if (errstat > 0) then
                call gdErrorHandler('EvaluateCostFunctionGSR: error encountered ' // &
                    'while evaluating goat equations. Setting cost function value ' // &
                    'to infinity and exiting evaluation', severityin=0)
                J = posinfval_R8()
                return
            end if

            ! Update goat (only if converged...)
            goat = goatproblem

        class default

            ! Throw error
            call gdErrorHandler('EvaluateCostFunctionGR: unexpected goat type')

        end select
        end associate



        ! Compute reduced cost function
        !==============================
        ! GOAT side
        call costfunction%costfunction%Evaluate(Jg, gradJg, hessJg, &
            goat, dogradient, dohessian, designvariables, &
            'goatvariables', goatvariables, gradJgoat)

        ! Associate coordinates for ease
        associate(x         => goat%grid%vert%x, &
                y         => goat%grid%vert%y)

        ! Evaluate magnetic field and derivatives
        call MFinterp%Evaluate(x, y, 0, 0, psi)
        call MFinterp%Evaluate(x, y, 1, 0, dpsidx)
        call MFinterp%Evaluate(x, y, 0, 1, dpsidy)

        ! Update SOLPS quantities
        costfunction%cfvsolps%g%vxX = x
        costfunction%cfvsolps%g%vxY = y
        costfunction%cfvsolps%g%vxFpsi = psi
        costfunction%cfvsolps%g%vxBx = dpsidx
        costfunction%cfvsolps%g%vxBy = dpsidy

        ! Update cell center coordinates (only necessary for vertex
        ! ordening, should be unnecessary in the future)
        do i = 1, goat%grid%cell%ntot
            tv = GetCellVert(goat%grid%cell, i)
            costfunction%cfvsolps%g%cvX(i) = sum(x(tv))/real(size(tv), kind=R8_G)
            costfunction%cfvsolps%g%cvY(i) = sum(y(tv))/real(size(tv), kind=R8_G)
        end do

        ! SOLPS side (gradient is w.r.t. coordinates, psi, dpsidx, dpsidy, ffbz)
        call costfunction%cfvsolps%Evaluate(Js, gradJs, dogradient)

        ! Total
        J = Jg + Js

        ! Compute goat linearization
        !===========================
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

        ! Compute gradients
        !==================
        ! Compute partial derivatives of psi and dpsidx, dpsidy w.r.t.
        ! grid coordinates (simply 'diagonal' matrix linearization
        ! (dpsidx, dpsidy computed before and up to date))
        call MFinterp%Evaluate(x, y, 2, 0, d2psidx2)
        call MFinterp%Evaluate(x, y, 1, 1, d2psidxdy)
        call MFinterp%Evaluate(x, y, 0, 2, d2psidy2)

        ! Compute solps gradient w.r.t. goat variables (assumed no contribution of ffbz)
        gradJsolps = gradJgoat ! initialize
        gradJsolps = 0
        gradJsolps(1:2*nv) = gradJs(1:2*nv) ! dJdx, dJdy
        gradJsolps(1:nv) = gradJsolps(1:nv) + gradJs(2*nv+1:3*nv)*dpsidx &
            + gradJs(3*nv+1:4*nv)*d2psidx2 + gradJs(4*nv+1:5*nv)*d2psidxdy ! dJdpsi*dpsidx + dJd(dpsidx)*d(dpsidx)dx + dJd(dpsidy)*d(dpsidy)dx
        gradJsolps(nv+1:2*nv) = gradJsolps(nv+1:2*nv) + gradJs(2*nv+1:3*nv)*dpsidy &
            + gradJs(3*nv+1:4*nv)*d2psidxdy + gradJs(4*nv+1:5*nv)*d2psidy2 ! dJdpsi*dpsidy + dJd(dpsidx)*d(dpsidx)dy + dJd(dpsidy)*d(dpsidy)dy

        ! Compute lagrange multipliers
        gradGgoat = jacGgoat%Transpose()
        call SolveSparseLinearSystemDI(gradGgoat, -(gradJgoat + gradJsolps), lambdaG, flag)

        ! Compute gradient
        gradGdes = jacGdes%Transpose()
        gradR = gradGdes%MatrixVectorProduct(lambdaG)
        gradJ = gradJg + gradR

        ! End association
        end associate

        ! Compute hessian
        !================
        if (dohessian) then
            ! Update the hessian approximation for the reduced part
            call costfunction%B%Update(designvariables%phi, gradR)

            ! Extract the hessian approximation and add with other
            ! contributions of which exact hessian is known
            hessJ = hessJg + costfunction%B%GetSparseHessian()
        end if

        ! Write out data for gradient verification
        !=========================================
        call WriteRealData('Js', Js)
        call WriteRealData('Jg', Jg)
        call WriteRealData('gradJ', gradJ)
        call WriteRealData('gradJg', gradJg)
        call WriteRealData('gradJs', gradJs)
        call WriteRealData('gradJsolps', gradJsolps)


        ! Housekeeping
        !=============
        ! Association
        end associate

        ! Optional arguments
        if (present(dJdvarin)) then
            dJdvarin = dJdvar
        end if
        if (present(dgradJdvarin)) then
            dgradJdvarin = dgradJdvar
        end if




    end subroutine


end module
