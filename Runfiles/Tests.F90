!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Short description
!==================
! This is the main driver to run tests for several components/features 
! of the goat module. Different tests are available 

module GOAT_tests 

    use mod_structured2Dgridding
    use, intrinsic :: ieee_arithmetic
    implicit none
    save 

    contains 

    !------------------------------------------------------------------!
    !                         General drivers                          !
    !------------------------------------------------------------------!
    ! Subroutine to run all tests
    subroutine RunAllTests()

        ! Interpolant testing
        !call TestStructuredInterpolant2D()
        !call TestCSparse()
        !call TestPLF2D
        !call TestQPSolvers
        !call TestDynamicArrays()
        !call TestContourTracing()
        !call TestSorting()
        call TestTopomeshGeneration()
    end subroutine

    !------------------------------------------------------------------!
    !                         Output routines                          !
    !------------------------------------------------------------------!
    ! Write starting header
    subroutine DisplayTestStart(testname)

        character(*),  intent(in)       :: testname 

        print *, '%---------------------------------------------------%'
        print *, '%                      TEST START                   %'
        print *, '%---------------------------------------------------%'
        print *, 'Test name: ' // testname 

    end subroutine

    ! Ending header
    subroutine DisplayTestEnd()

        print *, '%---------------------------------------------------%'
        print *, '%                      TEST END                     %'
        print *, '%---------------------------------------------------%'

    end subroutine

    !------------------------------------------------------------------!
    !                         Test routines                            !
    !------------------------------------------------------------------!
    ! Structured interpolant testing
    subroutine TestStructuredInterpolant2D()

        ! Modules
        !========
        use Interpolant
        use mod_sparseinterface

        implicit none
        save 

        ! Declare variables
        !==================
        ! Arguments

        ! Loop variables

        ! Auxiliary
        type(StructuredInterpolant2DUDT)    :: interp
        type(MySparseUDT)                   :: jacvq, jaca

        integer(I8)                         :: nx, ny, nq, nres, &
            derivx, derivy
        real(R8)                            :: Lx, Ly

        character(:), allocatable           :: meth 
        integer(I8), allocatable            :: reshuffle(:), aind(:), &
            vind(:)
        real(R8), allocatable               :: xgv(:), ygv(:), &
            v(:, :), a(:), vq(:, :), xq(:), yq(:), vqan(:, :), &
            dv(:, :), temp(:), relerr(:), vvals(:), avals(:), d(:), &
            ainit(:), vqeval(:), tempvals(:, :)

        real(R8), allocatable, dimension(:) :: vqFW, vqBW, dvqdaFW, &
            dvqdaBW, dvqdaC, eabsFW, eabsBW, eabsC, &
            erelFW, erelBW, erelC, dvqda

        ! Loop
        integer(I8)                         :: i, j, k 

        ! Initialize
        !===========
        ! Display
        call DisplayTestStart('TestStructuredInterpolant2D')

        ! Build structured grid
        nx = 10 ! number of cells, not vertices!
        ny = 10
        Lx = 10
        Ly = 5 
        allocate(xgv(nx+1), ygv(ny+1), v(nx+1, ny+1))
        xgv = Lx*[(k, k=0, nx)]/nx 
        ygv = Ly*[(k, k=0, ny)]/ny

        ! Set evaluation points (avoid out of bounds)
        !nq = 10000
        !nres = 6 ! number of results 
        !allocate(xq(nq), yq(nq), vqan(nq, nres), dv(nq, nres), &
        !    vq(nq, nres), temp(nq), vqeval(nq), vqFW(nq), vqBW(nq))
        !call random_number(xq)
        !call random_number(yq)

        nq = 2
        nres = 6 ! number of results 
        allocate(xq(nq), yq(nq), vqan(nq, nres), dv(nq, nres), vq(nq, nres), &
            temp(nq), vqeval(nq), vqFW(nq), vqBW(nq))
        xq = [0.12, 0.999]
        yq = [0.16, 0.999]
        
        xq = xq*Lx
        yq = yq*Ly

        ! Linear test case
        !=================
        ! Build interpolant values
        allocate(a(3)) ! coefficients
        a = [0, 1, 3]
        do i = 1, nx+1
            do j = 1, ny+1
                v(i, j) = a(1) + a(2)*xgv(i) + a(3)*ygv(j)
            end do 
        end do

        ! Build interpolant
        meth = 'uniformgrid'
        call interp%SetParameters(meth, 3, 6)
        call interp%Construct(xgv, ygv, v)


        ! Test evaluation
        !----------------
        ! Test1: evaluate the interpolant at some query points and 
        ! evaluate the analytic solution
        vqan(:, 1) = a(1) + a(2)*xq + a(3)*yq ! field value
        vqan(:, 2) = a(3) ! dfdy
        vqan(:, 3) = 0 ! d2fdy2
        vqan(:, 4) = a(2) !dfdx
        vqan(:, 5:6) = 0 ! d2fdx2, d2fdxdy

        k = 1
        do i = 0, 2
            do j = 0, 2-i
                call interp%Evaluate(xq, yq, i, j, temp)
                vq(:, k) = temp
                k = k + 1
            end do 
        end do 

        ! Evaluate norm
        dv = vqan - vq
        allocate(relerr(nres), reshuffle(nres))
        do k = 1, nres
            relerr(k) = maxval(abs(dv(:, k)/vqan(:, k)))
        enddo 
        reshuffle = [1, 4, 2, 6, 5, 3]

        ! Print
        print *, 'test case     f     dfdx    dfdy    d2fdx2      d2fdxdy     d2fdy2      '
        print *, 'Linear field test case', relerr(reshuffle)

        ! Test derivatives
        !-----------------
        ! Set fixed evaluation points to know in which cell we need to check
        ! (otherwise extremely likely to have zero derivatives everywhere)
        deallocate(xq, yq, vq, vqan, dv, temp, vqeval, vqFW, vqBW)
        nq = 2
        nres = 6 ! number of results 
        allocate(xq(nq), yq(nq), vqan(nq, nres), dv(nq, nres), vq(nq, nres), &
            temp(nq), vqeval(nq), vqFW(nq), vqBW(nq))
        xq = [0.12, 0.999]
        yq = [0.16, 0.999]
        
        xq = xq*Lx
        yq = yq*Ly

        ! Set derivatives
        derivx = 0
        derivy = 0

        ! Print
        print *, 'Partial derivatives, derivtype: values w.r.t. interpolant coefficients'
        
        ! Derivatives w.r.t. coefficients and initial values. Compare to FD
        call interp%EvaluateDiffInterp2Coef(xq, yq, derivx, derivy, vqeval, jaca)

        ! Compare with FD
        avals = reshape(interp%a, [size(interp%a, 1)*size(interp%a, 2)])

        ! Indices to test all coefficients in certain cell
        !aind = [(k, k = 12, 12+(size(interp%a, 2)-1)*nx*ny, nx*ny)] ![(k, k = 1, (nx+1)*(ny+1))]
        
        ! Indices to test only a couple of derivatives
        aind = [12, 12 + nx*ny, 12+2*nx*ny]

        d = [1e-8, 1e-6, 1e-4, 1e-2] 
        
        do i = 1, size(aind, 1)
            ! Output
            print *, 'variable: ', aind(i)
            print *, '| step size | eabsFW | eabsBW | eabsC | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call jaca%ExtractColumnFull(dvqda, aind(i))

            ! Loop over fd steps
            do j = 1, size(d, 1) 
                ! Forward difference
                avals(aind(i)) = avals(aind(i)) + d(j)
                interp%a = reshape(avals, [size(interp%a, 1), size(interp%a, 2)])
                call interp%Evaluate(xq, yq, derivx, derivy, vqFW)
                dvqdaFW = (vqFW - vqeval)/d(j)

                ! Backward difference
                avals(aind(i)) = avals(aind(i)) - 2*d(j)
                interp%a = reshape(avals, [size(interp%a, 1), size(interp%a, 2)])
                call interp%Evaluate(xq, yq, derivx, derivy, vqBW)
                dvqdaBW = (vqBW - vqeval)/(-d(j))

                ! Central difference
                dvqdaC = 0.5*(dvqdaFW + dvqdaBW)

                ! Reset value
                avals(aind(i)) = avals(aind(i)) + d(j)
                interp%a = reshape(avals, [size(interp%a, 1), size(interp%a, 2)])

                ! Errors
                eabsFW = abs(dvqda - dvqdaFW)
                eabsBW = abs(dvqda - dvqdaBW)
                eabsC = abs(dvqda - 0.5*(dvqdaFW + dvqdaBW))
                erelFW = eabsFW/dvqda 
                erelBW = eabsBW/dvqda 
                erelC = eabsC/dvqda 

                ! Print out information
                print *, d(j), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(abs(erelFW)), maxval(abs(erelBW)), &
                    maxval(abs(erelC)), maxloc(abs(erelC))

            end do
        end do

        ! Test derivatives
        !-----------------
        ! Print
        print *, 'Partial derivatives, derivtype: values w.r.t. interpolant coefficients'
        
        ! Derivatives w.r.t. coefficients and initial values. Compare to FD
        call interp%EvaluateDiffCoef2Val(xgv, ygv, v, jacvq)

        ! Compare with FD
        vvals = reshape(v, [(nx+1)*(ny+1)]) 
        ainit = reshape(interp%a, [size(interp%a)])
        vind = [1, nx + 1, 2, 2+ nx]
        d = [1e-8, 1e-6, 1e-4, 1e-2] 
        
        do i = 1, size(vind, 1)
            ! Output
            print *, 'variable: ', vind(i)
            print *, '| step size | eabsFW | eabsBW | eabsC | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call jacvq%ExtractColumnFull(dvqda, vind(i))

            ! Loop over fd steps
            do j = 1, size(d, 1) 
                ! Forward difference
                vvals(vind(i)) = vvals(vind(i)) + d(j)
                tempvals = reshape(vvals, [nx+1, ny+1])
                call interp%Construct(xgv, ygv, tempvals)
                dvqdaFW = (reshape(interp%a, [size(interp%a)]) - ainit)/d(j)

                ! Backward difference
                vvals(vind(i)) = vvals(vind(i)) - 2*d(j)
                tempvals = reshape(vvals, [nx+1, ny+1])
                call interp%Construct(xgv, ygv, tempvals)
                dvqdaBW = (reshape(interp%a, [size(interp%a)]) - ainit)/(-d(j))

                ! Central difference
                dvqdaC = 0.5*(dvqdaFW + dvqdaBW)

                ! Reset value
                vvals(vind(i)) = vvals(vind(i)) + d(j)

                ! Errors
                eabsFW = abs(dvqda - dvqdaFW)
                eabsBW = abs(dvqda - dvqdaBW)
                eabsC = abs(dvqda - 0.5*(dvqdaFW + dvqdaBW))
                erelFW = eabsFW/dvqda 
                erelBW = eabsBW/dvqda 
                erelC = eabsC/dvqda 

                ! Print out information
                print *, d(j), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(abs(erelFW)), maxval(abs(erelBW)), &
                    maxval(abs(erelC)), maxloc(abs(erelC))


            end do
        end do
        
        ! Test derivatives
        !-----------------
        ! Print
        print *, 'Partial derivatives, derivtype: values w.r.t. initial values coefficients'
        
        ! Set derivatives
        derivx = 0
        derivy = 0

        ! Compute values
        do i = 1, nx+1
            do j = 1, ny+1
                v(i, j) = sin(xgv(i))*sin(ygv(j))
            end do 
        end do

        ! Derivatives w.r.t. coefficients and initial values. Compare to FD
        call interp%Construct(xgv, ygv, v)
        call interp%Evaluate(xq, yq, derivx, derivy, vqeval)
        call interp%EvaluateDiffInterp2Val(xq, yq, derivx, derivy, jacvq)

        ! Compare with FD
        vvals = reshape(v, [(nx+1)*(ny+1)]) 
        vind = [1, nx + 1, 2, 2+ nx]
        d = [1e-8, 1e-6, 1e-4, 1e-2] 
        
        do i = 1, size(vind, 1)
            ! Output
            print *, 'variable: ', vind(i)
            print *, '| step size | eabsFW | eabsBW | eabsC | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call jacvq%ExtractColumnFull(dvqda, vind(i))

            ! Loop over fd steps
            do j = 1, size(d, 1) 
                ! Forward difference
                vvals(vind(i)) = vvals(vind(i)) + d(j)
                tempvals = reshape(vvals, [nx+1, ny+1])
                call interp%Construct(xgv, ygv, tempvals)
                call interp%Evaluate(xq, yq, derivx, derivy, vqFW)
                dvqdaFW = (vqFW - vqeval)/(d(j))

                ! Backward difference
                vvals(vind(i)) = vvals(vind(i)) - 2*d(j)
                tempvals = reshape(vvals, [nx+1, ny+1])
                call interp%Construct(xgv, ygv, tempvals)
                call interp%Evaluate(xq, yq, derivx, derivy, vqBW)
                dvqdaBW = (vqBW - vqeval)/(-d(j))

                ! Central difference
                dvqdaC = 0.5*(dvqdaFW + dvqdaBW)

                ! Reset value
                vvals(vind(i)) = vvals(vind(i)) + d(j)

                ! Errors
                eabsFW = abs(dvqda - dvqdaFW)
                eabsBW = abs(dvqda - dvqdaBW)
                eabsC = abs(dvqda - 0.5*(dvqdaFW + dvqdaBW))
                erelFW = eabsFW/dvqda 
                erelBW = eabsBW/dvqda 
                erelC = eabsC/dvqda 

                ! Print out information
                print *, d(j), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(abs(erelFW)), maxval(abs(erelBW)), &
                    maxval(abs(erelC)), maxloc(abs(erelC))

            end do
        end do

        ! Housekeeping
        deallocate(a)
        call interp%deallocate()

        ! Quadratic case
        !===============
        ! Build interpolant values
        allocate(a(6)) ! coefficients
        a = [0.0, 1.0, 3.0, 0.01, 0.05, 0.0016]
        do i = 1, nx+1
            do j = 1, ny+1
                v(i, j) = a(1) + a(2)*xgv(i) + a(3)*ygv(j) + &
                    a(4)*xgv(i)*ygv(j) + a(5)*xgv(i)**2 + a(6)*ygv(j)**2
            end do 
        end do

        ! Build interpolant
        meth = 'uniformgrid'
        call interp%SetParameters(meth, 3, 6)
        call interp%Construct(xgv, ygv, v)

        ! Test1: evaluate the interpolant at some query points and 
        ! evaluate the analytic solution
        vqan(:, 1) = a(1) + a(2)*xq + a(3)*yq + a(4)*xq*yq + &
            a(5)*xq**2 + a(6)*yq**2 ! field value
        vqan(:, 2) = a(3) + a(4)*xq + 2*a(6)*yq ! dfdy
        vqan(:, 3) = 2*a(6) ! d2fdy2
        vqan(:, 4) = a(2) + a(4)*yq + 2*a(5)*xq !dfdx
        vqan(:, 5) = a(4) !d2fdxdy
        vqan(:, 6) = 2*a(5) ! d2fdx2, d2fdxdy

        k = 1
        do i = 0, 2
            do j = 0, 2-i
                call interp%Evaluate(xq, yq, i, j, temp)
                vq(:, k) = temp
                k = k + 1
            end do 
        end do 

        ! Evaluate norm
        dv = vqan - vq
        do k = 1, nres
            relerr(k) = maxval(abs(dv(:, k)/vqan(:, k)))
        enddo 

        ! Print
        print *, 'Quadratic field test case', relerr(reshuffle)

        ! Housekeeping
        deallocate(a)
        call interp%deallocate()

        ! Housekeeping
        !=============
        deallocate(xgv, ygv, xq, yq, vqan, dv, temp, v)
        call DisplayTestEnd()



    end subroutine 

    ! CSparse interface
    subroutine TestCSparse()

        ! Description
        !============
        ! Test some CSparse implementation
        use Clayer
        use mod_sparseinterface

        ! Declare variables
        !==================
        ! Auxiliary
        type(CSparseUDT)   :: A2, A3
        type(MySparseUDT)  :: A, B, T, C, D, E, F
        integer(I8)        :: n
        real(c_double), allocatable     :: val(:)
        integer(c_int), allocatable     :: row(:), col(:)
        integer(c_int)                  :: nrow, ncol, nval 
        real(c_double), pointer         :: valp(:)
        integer(c_int), pointer         :: rowp(:), colp(:), rowp2(:)
        real(R8)                        :: t_start, t_end

        ! Create a sparse matrix
        nrow = 8
        ncol = 8
        nval = 5
        allocate(val(nval), row(nval), col(nval))
        val = [0.23, 1.23, 2.5, 4.0, 5.0]
        row = [1, 2, 3, 4, 5]
        col = [1, 2, 3, 4, 5]

        A%nrow = nrow 
        A%ncol = ncol 
        A%nval = nval 
        call A%Allocate()
        A%row = row 
        A%col = col 
        A%val = val 

        B = A*A 

        print *, B%row 
        print *, B%col 
        print *, B%val 

        ! Random sparse matrices
        !=======================
        ! Create
        n = 200000
        call CreateRandomSparseMatrix(D, 10000, 5000, n)
        call CreateRandomSparseMatrix(E, 5000, 10000, n)

        ! A = SpMMFBasic(D%nrow, D%ncol, D%nval, E%nrow, E%ncol, E%nval, D%row, D%col, D%val, E%row, E%col, E%val)
        !call c_f_pointer(A%col, rowp, [A%nval])
        !print *, rowp

        !print *, D%row 
        !print *, D%col 
        !print *, D%val

        ! Time and multiply
        call cpu_time(t_start)
        F = D*E 
        call cpu_time(t_end)
        print *, 'Time elapsed for multiplying matrices with ',  n, ' nonzeros: ', t_end-t_start 

    end subroutine

    ! Polygon Levelset functions
    subroutine TestPLF2D()

        ! Description
        !============
        ! Test whether the polygon level set routines are correctly set
        ! up. Here, we just write out some plot data to be compared...

        ! The usual
        use PolygonLevelsetFunction2D
        use mod_polygon
        use mod_sparseinterface

        ! Declare variables
        !==================
        ! Arguments

        ! Auxiliary
        class(PolygonLevelsetFunction2DUDT), allocatable    :: plfg, &
            plfce, plfca

        type(PLF2DGeneralOptionsUDT)                :: optionsg 
        type(PLF2DClosedExactOptionsUDT)            :: optionsce
        type(PLF2DClosedApproximationOptionsUDT)    :: optionsca 

        type(PolygonSetUDT)         :: psg, psce, psca

        type(MySparseUDT)           :: dplfgdvar, dplfcedvar, dplfcadvar

        integer(I8)                 :: nvaluesg, nxm, nym, derivx, derivy
        integer(I8), allocatable    :: varind(:)

        real(R8)                            :: NaN, maxxm, minxm, maxym, &
            minym
        real(R8), allocatable, dimension(:) :: xg, yg, xcps, &
            ycps, xncp, yncp, d, xp, yp, valuesg, vqg, xq, yq, xgv, ygv, &
            vqgFW, vqgBW, dvqFW, dvqBW, dvqC, eabsFW, eabsBW, eabsC, &
            erelFW, erelBW, erelC, dvqg, eabsdir, ereldir

        ! Loop
        integer(I8)                         :: i, j

        ! Initialize
        !===========
        ! Set NaN
        NaN = ieee_value(NaN, ieee_quiet_nan)

        ! Coordinates
        xg = [1, 1, 2, 4, -1]
        yg = [0, 2, 2, 1, 0]

        xcps = [1, 1, 2, 4, 1]
        ycps = [0, 2, 2, 1, 0]

        xncp = [1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.75, 0.75, 0.25, 0.25, 0.75 ]
        yncp = [0.0 , 1.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.75, 0.75, 0.25, 0.25 ]

        xncp(6) = NaN 
        yncp(6) = NaN

        ! Construct
        !==========
        ! Construct general polygon set
        call psg%Construct(xg, yg)

        ! Construct simple closed polygon set
        !call psce%Construct(xcps, ycps)
        call psce%Construct(xncp, yncp)

        ! Construct nested closed polygon set
        call psca%Construct(xncp, yncp)
        !call psca%Construct(xcps, ycps)

        ! Set plf options (only for ca)
        optionsca%C         = 3
        optionsca%M         = 6
        optionsca%meth      = 'uniformgrid'
        optionsca%offsetx   = 1
        optionsca%offsety   = 1
        optionsca%resx      = 20
        optionsca%resy      = 10
        optionsca%optionsClosedExact = optionsce

        ! Construct PLFs
        call InitializePolygonLevelsetFunction2D(plfg, psg, optionsg)
        call InitializePolygonLevelsetFunction2D(plfce, psce, optionsce)
        call InitializePolygonLevelsetFunction2D(plfca, psca, optionsca)

        ! Visualize
        call plfg%Visualize('generalplf')
        call plfce%Visualize('closedexactplf')
        call plfca%Visualize('closedapproximationplf')

        ! Test derivatives w.r.t. polygonset coordinates
        !===============================================
        ! General plf
        !------------
        ! Set finite differences & values
        d = [1e-8, 1e-6, 1e-4, 1e-2]
        varind = [1, 6, 2, 7, 3, 8, 4, 9, 5, 10]

        ! Extract polygonset coordinates
        call plfg%ps%GetVertices(xp, yp)
        nvaluesg = size(xp, 1)
        valuesg = [xp, yp]

        ! Create a meshgrid for polygon evaluation, limited by polygon
        ! extent
        nxm = 2
        nym = 2
        allocate(xq(nxm*nym), yq(nxm*nym), vqg(nxm*nym), vqgFW(nxm*nym), &
            vqgBW(nxm*nym))
        maxxm = maxval(xp)
        minxm = minval(xp)
        maxym = maxval(yp)
        minym = minval(yp)
        
        xgv = [(i, i = 0, nxm-1)]*(1./real(nxm-1, kind=R8))*(maxxm - minxm) + minxm 
        ygv = [(i, i = 0, nym-1)]*(1./real(nym-1, kind=R8))*(maxym - minym) + minym 

        call Construct2DStructuredGrid(xgv, ygv, nxm, nym, xq, yq)

        ! Evaluate implemented derivatives
        call plfg%Evaluate(xq, yq, &
            0, 0, vqg, 'polygonsetcoordinates', valuesg, dplfgdvar)
        
        ! Compute finite differences (crude)
        print *, 'Evaluating FD for general polygonset (variable: polygonsetcoordinates)'
        
        do j = 1, size(varind)
            ! Output
            print *, 'variable: ', varind(j)
            print *, '| step size | eabsdir | eabsFW | eabsBW | eabsC | ereldir | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call dplfgdvar%ExtractColumnFull(dvqg, varind(j))
            do i = 1, size(d)
                ! Forward
                !--------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i)

                ! Update polygonset coordinates
                call plfg%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfg%Initialize(plfg%ps)

                ! Compute
                call plfg%Evaluate(xq, yq, 0, 0, vqgFW)
                dvqFW = (vqgFW - vqg)/(d(i))

                ! Backward
                !---------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) - 2*d(i) ! 2*d to compensate for previous + d

                ! Update polygonset coordinates
                call plfg%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfg%Initialize(plfg%ps)

                ! Compute
                call plfg%Evaluate(xq, yq, 0, 0, vqgBW)
                dvqBW = (vqgBW - vqg)/(-d(i))

                ! Compute errors
                !---------------
                eabsFW = abs(dvqg - dvqFW)
                eabsBW = abs(dvqg - dvqBW)
                eabsC = abs(dvqG - 0.5*(dvqFW + dvqBW))
                erelFW = eabsFW/dvqg 
                erelBW = eabsBW/dvqg 
                erelC = eabsC/dvqg 
                eabsdir = eabsFW 
                ereldir = eabsFW/dvqg
                where (eabsFW > eabsBW)
                    eabsdir = eabsBW 
                    ereldir = eabsBW/dvqg
                end where
            

                ! Print out information
                print *, d(i), maxval(eabsdir), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(ereldir), maxval(erelFW), maxval(erelBW), &
                    maxval(erelC), maxloc(erelC)

                ! Update
                !-------
                ! Downdate values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i) ! + d to compensate for - d
            
                ! Update polygonset coordinates (no need to update plf)
                call plfg%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

            end do
        end do

        ! Housekeeping
        deallocate(xq, yq, vqg, vqgFW, vqgBW)

        

        ! Closed exact
        !-------------
        ! Set derivatives
        derivx = 0
        derivy = 0

        ! Set finite differences & values
        d = [1e-8, 1e-6, 1e-4, 1e-2]
        varind = [1, 5, 2, 6, 3, 7, 4, 8]

        ! Extract polygonset coordinates
        call plfce%ps%GetVertices(xp, yp)
        nvaluesg = size(xp, 1)
        valuesg = [xp, yp]

        ! Create a meshgrid for polygon evaluation, limited by polygon
        ! extent
        nxm = 33
        nym = 33
        allocate(xq(nxm*nym), yq(nxm*nym), vqg(nxm*nym), vqgFW(nxm*nym), &
            vqgBW(nxm*nym))
        maxxm = maxval(xp)+0.1
        minxm = minval(xp)-0.1
        maxym = maxval(yp)+0.1
        minym = minval(yp)-0.1
        
        xgv = [(i, i = 0, nxm-1)]*(1./real(nxm-1, kind=R8))*(maxxm - minxm) + minxm 
        ygv = [(i, i = 0, nym-1)]*(1./real(nym-1, kind=R8))*(maxym - minym) + minym 

        call Construct2DStructuredGrid(xgv, ygv, nxm, nym, xq, yq)

        ! Evaluate implemented derivatives
        call plfce%Evaluate(xq, yq, &
            derivx, derivy, vqg, 'polygonsetcoordinates', valuesg, dplfcedvar)

        ! Compute finite differences (crude)
        print *, 'Evaluating FD for closed exact polygonset (variable: polygonsetcoordinates)'
        
        do j = 1, size(varind)
            ! Output
            print *, 'variable: ', varind(j)
            print *, '| step size | eabsdir | eabsFW | eabsBW | eabsC | ereldir | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call dplfcedvar%ExtractColumnFull(dvqg, varind(j))
            do i = 1, size(d)
                ! Forward
                !--------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i)

                ! Update polygonset coordinates
                call plfce%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfce%Initialize(plfce%ps)

                ! Compute
                call plfce%Evaluate(xq, yq, derivx, derivy, vqgFW)
                dvqFW = (vqgFW - vqg)/d(i)

                ! Backward
                !---------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) - 2*d(i) ! 2*d to compensate for previous + d

                ! Update polygonset coordinates
                call plfce%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfce%Initialize(plfce%ps)

                ! Compute
                call plfce%Evaluate(xq, yq, derivx, derivy, vqgBW)
                dvqBW = (vqgBW - vqg)/(-d(i))

                ! Compute errors
                !---------------
                eabsFW = abs(dvqg - dvqFW)
                eabsBW = abs(dvqg - dvqBW)
                eabsC = abs(dvqG - 0.5*(dvqFW + dvqBW))
                eabsdir = eabsFW 
                ereldir = eabsFW/dvqg
                where (eabsFW > eabsBW)
                    eabsdir = eabsBW 
                    ereldir = eabsBW/dvqg
                end where
                
                erelFW = eabsFW/dvqg 
                erelBW = eabsBW/dvqg 
                erelC = eabsC/dvqg 
            

                ! Print out information
                print *, d(i), maxval(eabsdir), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(abs(ereldir)), maxval(abs(erelFW)), &
                    maxval(abs(erelBW)), maxval(abs(erelC)), maxloc(abs(erelC))


                ! Update
                !-------
                ! Downdate values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i) ! + d to compensate for - d
            
                ! Update polygonset coordinates (no need to update plf)
                call plfce%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

            end do
        end do

        ! Housekeeping
        deallocate(xq, yq, vqg, vqgFW, vqgBW, xp, yp)

        ! Closed approximation
        !---------------------
        ! Set derivatives
        derivx = 0
        derivy = 0

        ! Set finite differences & values
        d = [1e-8, 1e-6, 1e-4, 1e-2]
        !varind = [1, 5, 2, 6, 3, 7, 4, 8]
        varind = [1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16]

        ! Extract polygonset coordinates
        call plfca%ps%GetVertices(xp, yp)
        nvaluesg = size(xp, 1)
        valuesg = [xp, yp]

        ! Create a meshgrid for polygon evaluation, limited by polygon
        ! extent
        nxm = 6
        nym = 6
        allocate(xq(nxm*nym), yq(nxm*nym), vqg(nxm*nym), vqgFW(nxm*nym), &
            vqgBW(nxm*nym))
        maxxm = maxval(xp)+0.1
        minxm = minval(xp)-0.1
        maxym = maxval(yp)+0.1
        minym = minval(yp)-0.1
        
        xgv = [(i, i = 0, nxm-1)]*(1./real(nxm-1, kind=R8))*(maxxm - minxm) + minxm 
        ygv = [(i, i = 0, nym-1)]*(1./real(nym-1, kind=R8))*(maxym - minym) + minym 

        call Construct2DStructuredGrid(xgv, ygv, nxm, nym, xq, yq)

        ! Evaluate implemented derivatives
        call plfca%Evaluate(xq, yq, &
            derivx, derivy, vqg, 'polygonsetcoordinates', valuesg, dplfcadvar)

        ! Compute finite differences (crude)
        print *, 'Evaluating FD for closed approximation polygonset (variable: polygonsetcoordinates)'
        
        do j = 1, size(varind)
            ! Output
            print *, 'variable: ', varind(j)
            print *, '| step size | eabsdir | eabsFW | eabsBW | eabsC | ereldir | erelFW | erelBW | erelC | indrelC |'
            
            ! Extract value of implemented derivative
            call dplfcadvar%ExtractColumnFull(dvqg, varind(j))
            do i = 1, size(d)
                ! Forward
                !--------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i)

                ! Update polygonset coordinates
                call plfca%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfca%Initialize(plfca%ps)

                ! Compute
                call plfca%Evaluate(xq, yq, derivx, derivy, vqgFW)
                dvqFW = (vqgFW - vqg)/d(i)

                ! Backward
                !---------
                ! Update values
                valuesg(varind(j)) = valuesg(varind(j)) - 2*d(i) ! 2*d to compensate for previous + d

                ! Update polygonset coordinates
                call plfca%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

                ! Reconstruct the levelset function
                call plfca%Initialize(plfca%ps)

                ! Compute
                call plfca%Evaluate(xq, yq, derivx, derivy, vqgBW)
                dvqBW = (vqgBW - vqg)/(-d(i))

                ! Compute errors
                !---------------
                eabsFW = abs(dvqg - dvqFW)
                eabsBW = abs(dvqg - dvqBW)
                eabsC = abs(dvqG - 0.5*(dvqFW + dvqBW))
                erelFW = eabsFW/dvqg 
                erelBW = eabsBW/dvqg 
                erelC = eabsC/dvqg 
                eabsdir = eabsFW 
                ereldir = eabsFW/dvqg
            
                where (eabsFW > eabsBW)
                    eabsdir = eabsBW 
                    ereldir = eabsBW/dvqg
                end where

                ! Print out information
                print *, d(i), maxval(eabsdir), maxval(eabsFW), maxval(eabsBW), &
                    maxval(eabsC), maxval(ereldir), maxval(abs(erelFW)), maxval(abs(erelBW)), &
                    maxval(abs(erelC)), maxloc(abs(erelC))

                ! Update
                !-------
                ! Downdate values
                valuesg(varind(j)) = valuesg(varind(j)) + d(i) ! + d to compensate for - d
            
                ! Update polygonset coordinates (no need to update plf)
                call plfca%ps%UpdateCoordinates(valuesg(1:nvaluesg), valuesg(nvaluesg+1:2*nvaluesg))

            end do
        end do

        ! Housekeeping
        deallocate(xq, yq, vqg, vqgFW, vqgBW, xp, yp)



    end subroutine

    ! QP solvers
    subroutine TestQPSolvers()

        ! Description
        !============
        ! Test the different QP solvers by solving unconstrained, 
        ! equality constrained, and inequality constrained problems. 
        ! We take a very simple 2D quadratic problem (quadratic cost 
        ! function, linear (in)equality constraints) with known 
        ! solution. Normally, the equality constrained problems with 
        ! direct solver should yield the exact solution up to machine
        ! precision (we provide the exact hessian, which is assumed 
        ! by the QP - otherwise we should do SQP/use a different solver)
        !
        ! The test problem is: 
        !
        !   min_x1,x2   0.5*( (x1 - x1*)^2 + (x2 - x2*)^2 ) 
        !                   = 0.5* ( x1^2 + x2^2) - (x1x1* + x2x2*) + c
        !   s.t.        (x1 - x1*) + (x2 - x2*) = 1
        !               (x1 - x1*) - (x2 - x2*) <= -1
        !
        ! Note that we neglect the constant term c. 
        ! For the unconstrained problem, the optimum lies at x1*, x2*.
        ! For the equality constrained problem, lambda = -0.5 and 
        ! x1 = 0.5 + x1*, x2 = 0.5 + x2*
        ! For the full problem, x1 = x1*, x2 = 1 + x2*, lambda = -0.5, 
        ! mu = 0.5

        ! Modules & the usual
        use optmod_qp
        use mod_constants
        use mod_sparseinterface
        use mod_precision
        use optmod_hessianapproximation

        implicit none 

        ! Declare variables
        !==================
        ! Arguments

        ! Auxiliary
        integer(I8)                             :: flag, maxit, &
            verbosity

        real(R8)                                :: x1s, x2s, tol 
        real(R8), allocatable, dimension(:)     :: gradJ, b, c, x, &
            lambda, mu, x0, lambda0, mu0, xe, lambdae, mue
        real(R8), allocatable, dimension(:, :)  :: Bdinit, jacGd, jacHd

        type(MySparseUDT)                       :: Bspinit, jacGsp, &
            jacHsp
        class(HessianApproximationUDT), allocatable     :: Bd, Bsp

        ! Initialize
        !===========
        ! Display
        call DisplayTestStart('TestQPSolvers')

        ! Construct test problem
        !=======================
        ! Problem parameters
        x1s = 1.0
        x2s = 2.0

        gradJ = [-x1s, -x2s]
        b = [1 + x1s + x2s]
        c = [-1 + x1s - x2s]

        ! Solver parameters
        maxit = 10
        tol = 1e-8
        verbosity = 2

        ! Initial guess
        x0 = [0, 0]*1.0
        lambda0 = [0]*1.0
        mu0 = [0]*1.0

        ! Dense representation
        Bdinit = reshape([1, 0, 0, 1]*1.0, [2, 2])
        Bd = ConstructHessianApproximation('no', 2, Bdinit)  

        jacGd = reshape([1, 1]*1.0, [1, 2])
        jacHd = reshape([1, -1]*1.0, [1, 2])

        ! Sparse representation
        Bspinit = ConstructMySparse(Bdinit)
        Bsp = ConstructHessianApproximation('no', 2, Bspinit)

        jacGsp = ConstructMySparse(jacGd)
        jacHsp = ConstructMySparse(jacHd)

        ! Unconstrained problem
        !----------------------
        ! Print
        print *, ' '
        print *, 'testing unconstrained problem'

        ! Set analytical solution
        xe = [x1s, x2s]

        ! Solve dense
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bd, gradJ, x, flag)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (dense): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)

        ! Solve sparse
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bsp, gradJ, x, flag)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (sparse): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)

        ! Equality constrained
        !---------------------
        ! Print
        print *, ' '
        print *, 'testing equality constrained problem'
        

        ! Set analytical solution
        xe = [x1s + 0.5, x2s + 0.5]
        lambdae = [-0.5]

        ! Solve dense
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bd, gradJ, jacGd, b, x, lambda, flag)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (dense): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution lambda (dense): ', &
            maxval(abs(lambda - lambdae)), maxval(abs(lambda - lambdae)/lambdae)

        ! Solve sparse
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bsp, gradJ, jacGsp, b, x, lambda, flag)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (sparse): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution lambda (sparse): ', &
            maxval(abs(lambda - lambdae)), maxval(abs(lambda - lambdae)/lambdae)

        ! Inequality constrained
        !-----------------------
        ! Print
        print *, ' '
        print *, 'testing equality constrained problem'

        ! Set analytical solution
        xe = [x1s, x2s + 1]
        lambdae = [-0.5]
        mue = [0.5]

        ! Solve dense
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bd, gradJ, jacGd, b, jacHd, c, x, lambda, &
            mu, flag, maxit, tol, verbosity)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (dense): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution lambda (dense): ', &
            maxval(abs(lambda - lambdae)), maxval(abs(lambda - lambdae)/lambdae)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution mu (dense): ', &
            maxval(abs(mu - mue)), maxval(abs(mu - mue)/mue)

        ! Solve sparse
        x = x0
        lambda = lambda0 
        mu = mu0 
        call SolveQPDirect(Bsp, gradJ, jacGsp, b, jacHsp, c, x, lambda, &
            mu, flag, maxit, tol, verbosity)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution (sparse): ', maxval(abs(x - xe)), maxval(abs(x - xe)/xe)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution lambda (sparse): ', &
            maxval(abs(lambda - lambdae)), maxval(abs(lambda - lambdae)/lambdae)
        print *, 'max absolute and relative difference between ' // & 
            'analytical and numerical solution mu (dense): ', &
            maxval(abs(mu - mue)), maxval(abs(mu - mue)/mue)


        ! Housekeeping
        !=============
        ! Display
        call DisplayTestEnd()

    end subroutine

    ! Dynamic array
    subroutine TestDynamicArrays()

        ! Modules
        !========
        use mod_precision
        use mod_dynamicarrays

        ! Declare variables
        !==================
        ! Auxiliary
        type(RealDynamicArrayUDT)       :: rda1, rda2, rda3
        real(R8), allocatable           :: val(:)

        ! Initialize
        !===========
        ! Call header
        call DisplayTestStart('Dynamic array testing')

        ! Tests
        !======
        ! Construct arrays
        val = [1, 2, 3, 4, 5]*1.0_R8
        rda1 = ConstructRealDynamicArray(val)
        rda2 = ConstructRealDynamicArray(val+2.5)

        ! Print
        print *, 'rda1 array values: ', rda1%val 
        print *, 'rda2 array values: ', rda2%val 

        ! Sum arrays
        rda3 = rda1 + rda2
        print *, 'sum of rda1 and rda2: ', rda3%val
        rda1 = rda1 + 1.0_R8
        print *, 'sum of rda1 and scalar: ', rda1%val
        rda1 = rda1 + val
        print *, 'sum of rda1 and array: ', rda1%val

        ! Extend arrays
        call rda1%Insert(2.5_R8, 4_I8)
        print *, 'rda1 extended values: ', rda1%val
        call rda2%Insert([2.0_R8, 3.0_R8], int([2, 6], kind=I8))
        print *, 'rda2 extended values: ', rda2%val

        ! Append arrays
        call rda1%Append(1.0_R8)
        call rda1%Append([2.0, 3.0]*1.0_R8)
        print *, 'rda1 appended values: ', rda1%val

        ! Remove array values


        

        ! Finalize
        !=========
        call DisplayTestEnd()

    end subroutine

    ! Contour tracing
    subroutine TestSorting()

        ! Description
        !============
        ! Test contour tracing algorithm for simple contours. Data is 
        ! written out in polygon format, to be plotted using python.

        ! Modules
        !========
        use mod_precision 
        use mod_constants
        use mod_plotter
        use mod_sort
        
        ! Declare variables
        !==================
        ! Auxiliary
        integer(I8), allocatable        :: int_a_rng(:), temp(:), &
            ind(:)
        integer(I8)                     :: n 

        real(R8), allocatable           :: real_a_rng(:), real_a_smallval(:)

        ! Initialize
        !===========
        ! Construct test matrices
        n = 10 ! start small

        ! Construct random vector of integers
        allocate(real_a_rng(n), int_a_rng(n), ind(n), real_a_smallval(n))
        call random_number(real_a_rng) 
        int_a_rng = floor(real_a_rng*n)
        real_a_smallval = 1e-13

        ! Print 
        print *, 'Unsorted integer array: ', int_a_rng 
        print *, 'Unsorted real array: ', real_a_rng

        ! Sort
        !=====
        temp = int_a_rng 
        call Sort(int_a_rng, ind=ind)
        call Sort(real_a_rng)
        call Sort(real_a_smallval)

        ! Print
        !======
        print *, 'Sorted integer array: ', int_a_rng
        print *, 'Indices: ', ind 
        print *, 'Sorted array by indexing', temp(ind)
        print *, 'Sorted real array: ', real_a_rng
        print *, 'Sorted small value array: ', real_a_smallval

    end subroutine

    ! Sorting
    subroutine TestContourTracing()

        ! Description
        !============
        ! Test contour tracing algorithm for simple contours. Data is 
        ! written out in polygon format, to be plotted using python.

        ! Modules
        !========
        use mod_precision 
        use mod_contour2D
        use mod_constants
        use mod_plotter
        use mod_dynamicarrays
        
        ! Declare variables
        !==================
        ! Auxiliary
        real(R8)                        :: Lx, Ly 
        real(R8), allocatable           :: xgv(:), ygv(:), xg(:), yg(:), &
            vtest(:, :), cval(:), ea(:), xw(:), yw(:)
        integer(I8)                     :: nx, ny
        integer(I8), allocatable        :: order(:), eai(:)
        logical, allocatable            :: issaddle(:)
        type(ContourUDT), allocatable   :: contours(:)
        type(RealDynamicArrayUDT)       :: xc, yc 

        ! Loop
        integer(I8)                     :: k 

        ! Initialize
        !===========
        ! Set the grid coordinates
        Lx = 1
        Ly = 1
        nx = 101
        ny = 201 
        xgv = real([(k, k = 0, nx-1)], kind=R8)*(Lx/real((nx-1), kind=R8))
        ygv = real([(k, k = 0, ny-1)], kind=R8)*(Ly/real((ny-1), kind=R8))

        ! Construct 2D grid
        allocate(xg(nx*ny), yg(nx*ny))
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Construct test values
        vtest = reshape(sin(xg*pi_R8)*sin(yg*pi_R8) &
            + sin(2*xg*pi_R8)*sin(2*yg*pi_R8) &
            + sin(3*xg*pi_R8)*sin(3*yg*pi_R8), [nx, ny])

        ! Trace contours
        !===============
        ! Set contour values
        cval = [0.2, 0.8, -0.3]*1.0_R8

        ! Set saddle points (empty arrays)
        allocate(ea(0), order(0), eai(0))
        allocate(issaddle(0))

        ! Trace
        call TraceContoursStructured2D(vtest, xgv, ygv, cval, ea, ea, ea, &
            eai, issaddle, order, contours)

        ! Visualize
        !==========
        ! Write as polygons
        xc = ConstructRealDynamicArray(contours(1)%x)
        yc = ConstructRealDynamicArray(contours(1)%y)
        do k = 2, size(contours)
            call xc%Append([nanval_R8(), contours(k)%x])
            call yc%Append([nanval_R8(), contours(k)%y])
        end do 
        xw = xc%Get()
        yw = yc%Get()
        call Write2DPolygonData(xw, yw, 'testcontourtracing')


    end subroutine

    ! Topological mesh generation
    subroutine TestTopomeshGeneration()

        ! Description
        !============
        ! Test topological mesh generation algorithm for predefined
        ! magnetic field and vessel geometry. Only the true necessary
        ! fields of the magnetic field and vessel are initialized here

        ! Modules
        !========
        use mod_precision 
        use mod_contour2D
        use mod_constants
        use mod_plotter
        use mod_dynamicarrays
        use goatmod_types 
        use ggmod_topology2D 
        use goatmod_userinput
        use ggmod_gridgeneration2D
        use mod_streamlinetracing2D
        
        ! Declare variables
        !==================
        ! Auxiliary
        character(:), allocatable       :: meth 
        real(R8)                        :: Lx, Ly 
        real(R8), allocatable           :: xgv(:), ygv(:), xg(:), yg(:), &
            vtest(:, :), cval(:), ea(:), xw(:), yw(:), xp(:), yp(:), &
            Vfx(:), Vfy(:), Vv(:)
        integer(I8)                     :: nx, ny
        integer(I8), allocatable        :: order(:) 
        type(ContourUDT), allocatable   :: contours(:)
        class(ContourTracerUDT), allocatable    :: fieldtracer, vesseltracer
        class(StreamlineTracerUDT), allocatable :: streamlinetracer
        type(RealDynamicArrayUDT)       :: xc, yc 
        type(magneticFieldUDT)          :: magneticField 
        type(VesselUDT)                 :: vessel 
        type(TopomeshUDT)               :: topomesh 
        type(GridUDT)                   :: simgrid
        type(PLF2DClosedExactOptionsUDT)    :: plfoptions
        type(TopomeshOptionsUDT)        :: topoptions
        type(GGOptionsUDT)              :: ggoptions
        real(R8), parameter             :: emptyR8(0)= 0
        integer(I8), parameter          :: emptyI8(0) = 0
        logical, parameter              :: emptyL(0) = .false.

        ! Loop
        integer(I8)                     :: k 

        ! Initialize
        !===========
        ! Set the vessel polygon coordinates (closed polygon)
        xp = [0.25, 0.75, 0.75, 0.25, 0.25]
        yp = [0.25, 0.25, 0.75, 0.75, 0.25]

        ! Set the grid coordinates
        Lx = 1
        Ly = 1
        nx = 101
        ny = 201 
        xgv = real([(k, k = 0, nx-1)], kind=R8)*(Lx/real((nx-1), kind=R8))
        ygv = real([(k, k = 0, ny-1)], kind=R8)*(Ly/real((ny-1), kind=R8))

        ! Construct 2D grid
        allocate(xg(nx*ny), yg(nx*ny))
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Construct test values
        vtest = reshape(sin(xg*pi_R8)*sin(yg*pi_R8) &
            + sin(2*xg*pi_R8)*sin(2*yg*pi_R8) &
            + sin(3*xg*pi_R8)*sin(3*yg*pi_R8), [nx, ny])

        ! Construct magnetic field
        meth = 'uniformgrid'
        call magneticField%interp%SetParameters(meth, 3, 6)
        call magneticField%interp%Construct(xgv, ygv, vtest)
        allocate(Vfx(nx*ny), Vfy(nx*ny), Vv(nx*ny))
        call magneticField%interp%Evaluate(xg, yg, 1, 0, Vfx)
        call magneticField%interp%Evaluate(xg, yg, 0, 1, Vfy)

        ! Construct vessel polygon set (exact representation for now)
        call vessel%polygonset%Construct(xp, yp)
        call InitializePolygonLevelsetFunction2D(vessel%plfvessel, vessel%polygonset, &
            plfoptions)
        call vessel%plfvessel%Evaluate(xg, yg, 0, 0, Vv)

        ! Set the default topomesh options
        call topoptions%SetDefaults()
        topoptions%vresx = 400
        topoptions%vresy = 400

        ! Set the default grid generation options
        call ggoptions%SetDefaults()

        ! Construct tracers
        !==================
        fieldtracer = ConstructStructuredTracer(vtest, xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topoptions%npmin, &
            topoptions%npmax, topoptions%dl)
        vesseltracer = ConstructStructuredTracer(&
            reshape(Vv, [nx, ny]), xgv, ygv, &
            emptyR8, emptyR8, emptyR8, emptyI8, emptyL, topoptions%npmin, &
            topoptions%npmax, topoptions%dl)
        streamlinetracer = ConstructStructuredStreamlineTracer(&
            reshape(Vfx, [nx, ny]), reshape(Vfy, [nx, ny]), xgv, ygv, &
            step=ggoptions%orthtracerstep, nsteps=ggoptions%orthtracernsteps)

        ! Construct topological mesh
        !===========================
        call ConstructTopologicalMesh(vessel, magneticField, topoptions, &
            topomesh, fieldtracer, vesseltracer, streamlinetracer)      
            
        ! Construct grid
        !===============
        call GenerateUnstructuredAlignedGrid(simgrid, topomesh, magneticField, &
            vessel, fieldtracer, vesseltracer, streamlinetracer, ggoptions)

    end subroutine


    !------------------------------------------------------------------!
    !                           Auxiliary                              !
    !------------------------------------------------------------------!

    subroutine CreateRandomSparseMatrix(A, nrow, ncol, nnz)

        ! Description
        !============
        ! Create a random sparse matrix with number of rows equal to 
        ! nrow (analogously ncol) and a number of nonzeros nnz. All is 
        ! randomly generated using the built in function random_number.
        ! The values are distributed over the interval [-1, 1]. 

        ! Modules
        !========
        use mod_sparseinterface

        ! Declare
        !========
        ! Arguments
        type(MySparseUDT)       :: A 
        integer(I8), intent(in) :: nrow, ncol, nnz 

        ! Auxiliary
        real(R8)                :: temp(1:nnz), val(1:nnz)
        integer(I8)             :: row(1:nnz), col(1:nnz)

        ! Test
        !======
        if (allocated(A%val)) then 
            call A%Deallocate()
        end if

        ! Construct
        !==========
        ! Rows
        call random_number(temp)
        row = ceiling(temp*nrow)

        ! Columns
        call random_number(temp)
        col = ceiling(temp*ncol)

        ! Values
        call random_number(temp)
        val = temp*2-1 ! to have some zero and nonzero values

        ! Assign
        A%nrow = nrow 
        A%ncol = ncol
        A%nval = nnz 
        A%row = row 
        A%col = col 
        A%val = val 



    end subroutine


end module 




program Tests

    ! Description
    !============
    ! Tests for the structured 2D interpolant. We test for a given 
    ! analytical field whether it is correctly represented again by 
    ! the structured polynomial interpolant

    use GOAT_tests 

    call RunAllTests()


end program Tests


                     
