!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the (abstract) PolygonLevelsetFunction2D type 
! along with other (non-abstract) derived types. The main object 
! provides (at least) evaluation routines for (x, y) coordinate pairs. 
! Each levelset is based on a polygon set object, where the zero level 
! set corresponds to the (approximation of) the polygon set. Note that 
! most levelset functions here are only defined for closed, 
! non-intersecting (possibly nested) sets of polygons. 

! For mathematical descriptions, see the definition of each extended
! type. Basically, there are three types: 
!
!   1)  the 'general' one, which is valid for any polygon 
!       (self-intersecting, closed/open, ...) but which 
!       is discontinuous. The description is exact. 
!   2)  a smoother (but still non-smooth in several regions) 
!       approximation for nested closed non-intersecting polygons only. 
!       Still an exact representation of the original polygon set
!   3)  a fully smooth yet approximate representation of 2) by using 
!       an interpolation function (hence this module requires the 
!       Interpolant module)

module PolygonLevelsetFunction2D 

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_constants
    use mod_plotter
    use Interpolant
    use mod_polygon
    use mod_sparseinterface
    use mod_structured2Dgridding
    use omp_lib

    ! The usual
    implicit none
    save
    private 

    ! Types
    public PolygonLevelsetFunction2DUDT, PolygonLevelsetFunction2DGeneralUDT, &
        PolygonLevelsetFunction2DClosedExactUDT, &
        PolygonLevelsetFunction2DClosedApproximationUDT
    public PLF2DOptionsUDT, PLF2DGeneralOptionsUDT, &
        PLF2DClosedExactOptionsUDT, PLF2DClosedApproximationOptionsUDT

    ! Routines
    public InitializePolygonLevelsetFunction2D, InitializePLF2DGeneral, &
        InitializePLF2DClosedExact, InitializePLF2DClosedApproximation
    public assignment(=)

    !==================================================================!
    !                                                                  !
    !                              TYPES                               !
    !                                                                  !
    !==================================================================! 

    ! Abstract
    !=========
    ! Options type
    type, abstract :: PLF2DOptionsUDT 

        ! Description
        !============
        ! Abstract options type to get initialization routine to work (?)
        character(:), allocatable   :: levelsettype 

    end type

    ! Polygon levelset function type
    type, abstract :: PolygonLevelsetFunction2DUDT 

        ! Description
        !============
        ! Main abstract type

        ! Fields
        type(PolygonSetUDT)     :: ps 

    contains

        ! Initialization
        procedure(InitializePLF2DINT), deferred     :: Initialize

        ! Evaluation
        procedure(EvaluatePLF2DINT), deferred       :: Evaluate 

        ! Visualization
        procedure   :: Visualize        => VisualizePolygonLevelsetFunction2D

    end type

    ! Extended
    !=========
    ! Options for general levelset function
    type, extends(PLF2DOptionsUDT)  :: PLF2DGeneralOptionsUDT

    contains 

        ! Options setter
        procedure :: Set    => SetOptionsGeneral

    end type

    ! Options for closed exact representation
    type, extends(PLF2DOptionsUDT)  :: PLF2DClosedExactOptionsUDT
    contains 

        ! Options setter
        procedure :: Set    => SetOptionsClosedExact

    end type

    ! Options for closed approximate representation
    type, extends(PLF2DOptionsUDT)  :: PLF2DClosedApproximationOptionsUDT

        ! Fields
        ! Options for closed approximation representation 
        type(PLF2DClosedExactOptionsUDT) :: optionsClosedExact

        ! Interpolation options
        real(R8)                    :: offsetx, offsety 
        real(R8), allocatable       :: xrange(:), yrange(:) ! determines extent of interpolant
        integer(I8)                 :: resx, resy, C, M 
        character(:), allocatable   :: meth 

    contains 

        ! Options setter
        procedure :: Set    => SetOptionsClosedApproximation

    end type
    
    ! General levelset function
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DGeneralUDT 

        ! Description
        !============
        ! General levelset type. The following fields are added for 
        ! easier computation later on (computed based on given polygon
        ! set):
        !
        ! - nc:         number of polygon edgees
        ! - x, y        nc-by-2 array, polygon egde coordinates ( first
        !               column is first node, second column is second 
        !               node of the polygon) 
        ! - xf, yf      nc-by-1 array with the polygon edge centers
        ! - nxp, nyp    the normalized polygon edge normals
        ! - txp, typ    the NOT normalized polygon tangents

        ! Mathematical description
        !=========================
        ! The mathematical formulation is based on two key requirements for a 
        ! point to strictly lie on a polygon edge: 
        !
        ! - The perpendicular distance between the line formed by the edge 
        ! points and the query point should be zero
        ! - The lengths of the vectors between the edge points and the query 
        ! point should equal the length of the vector between the edge points
        ! 
        ! The function value can therefore be defined as (for one edge, one 
        ! sample point):
        !
        !       f = abs (abs(vxp1) + abs(vxp2) - abs(txp)) 
        !           + abs (abs(vyp1) + abs(vyp2) - abs(typ))
        !           + abs(dvn)
        ! 
        ! Here, vp1, vp2 are the vectors between the query point and the first
        ! and second edge points, tp is the vector between the two edge points
        ! (and their components are indicated with 'x' and 'y'). Dvn is the 
        ! perpendicular length between the query point and the line passing 
        ! through the edge points, which is computed by projecting the vector
        ! between the edge center (xf, yf) and query point with the edge normal:
        !
        !       dvn = vx * nxp + vy * nyp
        !
        ! where vx, vy are the components of the vertex - edge face vector. 

        ! For the case with multiple edges, we simply loop over all edges and
        ! compute the minimal value of this function for each of the edges. 
        ! Though the abs function only results in non-smoothness at the edge,
        ! the max function leads to non-smoothness at locations away from the 
        ! edge.  

        ! Fields
        type(PLF2DGeneralOptionsUDT)    :: options
        integer(I8)                     :: nc
        real(R8), allocatable           :: xf(:), yf(:), nxp(:), &
            nyp(:), txp(:), typ(:), x(:, :), y(:, :)

    contains 

        ! Initializer
        procedure :: Initialize     => InitializePLF2DGeneral

        ! Evaluate
        procedure :: Evaluate       => EvaluatePLF2DGeneral

    end type

    ! Closed exact representation
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DClosedExactUDT 

        ! Description
        !============
        ! Exact representation of closed, non-intersecting set of 
        ! polygons. Hereto, we subdivide the 2D domain into vertex and
        ! edge regions, which then determines how the distance to the 
        ! polygon is computed (either the Euclidean distance to the 
        ! vertex, or the normal distance to the line formed by a 
        ! certain polygon edge). 

        ! Note: an extra method is provided specifically for this 
        ! routine to determine from which edge and which polygon 
        ! the value was interpolated. Useful for some reconstruction
        ! purposes (e.g. determining face labels based on vessel 
        ! structures). It should be noted that at polygon transitions, 
        ! there will be a discontinuity in value (on the vertex itself, 
        ! two values are actually present), and we simply take one 
        ! value in the vertex itself. 

        ! Fields
        type(PLF2DClosedExactOptionsUDT)    :: options
        real(R8), allocatable               :: xp(:), yp(:), xf(:), &
            yf(:), nxp(:), nyp(:), tnp(:), nxpv(:, :), &
            nypv(:, :), crossprod(:), thetav(:), xp1(:), xp2(:), &
            yp1(:), yp2(:)
        integer(I8), allocatable, dimension(:)      :: vp, vp1, vp2
        integer(I8), allocatable, dimension(:, :)   :: edgelabel, &
            vertlabel

    contains 

        ! Initialization
        procedure :: Initialize     => InitializePLF2DClosedExact

        ! Evaluate
        procedure :: Evaluate       => EvaluatePLF2DClosedExact
        procedure :: EvaluateLabel  => EvaluatePLF2DClosedExactLabel
        
        ! Visualization of label 
        procedure :: VisualizeLabel     => VisualizePolygonLevelsetFunctionLabelClosedExact

    end type

    ! Closed approximate representation
    type, extends(PolygonLevelsetFunction2DUDT)  :: &
        PolygonLevelsetFunction2DClosedApproximationUDT 

        ! Options
        type(PLF2DClosedApproximationOptionsUDT)    :: options

        ! Interpolant
        type(StructuredInterpolant2DUDT)            :: interp

    contains 

        ! Initialization
        procedure :: Initialize     => InitializePLF2DClosedApproximation

        ! Evaluation
        procedure :: Evaluate       => EvaluatePLF2DClosedApproximation

        ! Grid construction
        procedure :: ConstructGrid  => ConstructGridPLF2DClosedApproximation

    end type

    
    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================! 

    ! Abstract
    !=========
    abstract interface 

        ! Initialization
        subroutine InitializePLF2DINT(plf, ps)
            import :: PolygonLevelsetFunction2DUDT, PolygonSetUDT, &
                PLF2DOptionsUDT
            class(PolygonLevelsetFunction2DUDT) :: plf 
            type(PolygonSetUDT), intent(in)     :: ps
        end subroutine

        ! Evaluation
        subroutine EvaluatePLF2DINT(plf, xq, yq, derivx, derivy, vq, &
            varin, valuesin, dplfdvarin)
            import :: PolygonLevelsetFunction2DUDT, R8, I8, MySparseUDT
            class(PolygonLevelsetFunction2DUDT)     :: plf 
            real(R8), intent(in)                    :: xq(:), yq(:)
            real(R8), intent(out)                   :: vq(size(xq))
            integer(I8), intent(in)                 :: derivx, derivy

            character(*), intent(in), optional      :: varin 
            real(R8), intent(in), optional          :: valuesin(:)
            type(MySparseUDT), optional             :: dplfdvarin


        end subroutine

    end interface

    ! Operator overloading
    interface assignment(=)
        module procedure AssignPLF2DClass
    end interface
    
    contains

    !==================================================================!
    !                                                                  !
    !                            ROUTINES                              !
    !                                                                  !
    !==================================================================! 

    !------------------------------------------------------------------!
    !                       ABSTRACT TYPE ROUTINES                     !
    !------------------------------------------------------------------!

    ! Initializer
    recursive subroutine InitializePolygonLevelsetFunction2D(plf, ps, options) 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DUDT), allocatable, intent(out)  :: plf 
        type(PolygonSetUDT), intent(in)                         :: ps
        class(PLF2DOptionsUDT), intent(in)                      :: options

        ! Initialize
        !===========
        if (allocated(plf)) then 
            deallocate(plf)
        end if 
        select type (options)

        type is (PLF2DGeneralOptionsUDT)

            ! Allocate
            allocate(PolygonLevelsetFunction2DGeneralUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DGeneralUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        type is (PLF2DClosedExactOptionsUDT) 

            ! Allocate
            allocate(PolygonLevelsetFunction2DClosedExactUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DClosedExactUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        type is (PLF2DClosedApproximationOptionsUDT)

            ! Allocate
            allocate(PolygonLevelsetFunction2DClosedApproximationUDT::plf)

            ! Initialize 
            select type(plf)

            type is (PolygonLevelsetFunction2DClosedApproximationUDT)

                plf%options = options
                call plf%Initialize(ps)

            end select

        end select 

    end subroutine

    ! Visualization
    subroutine VisualizePolygonLevelsetFunction2D(plf, savefilepath, &
            nxin, nyin)

        ! Description
        !============
        ! This routine provides a crude way of visualizing the 
        ! levelset function by making a 2D patchplot on a structured 
        ! grid. The extent is equal to the coordinate boundaries 
        ! of the polygon set vertices, and the resolution in both 
        ! directions is fixed and constant. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DUDT)      :: plf 
        character(*), intent(in)                 :: savefilepath
        integer(I8), intent(in), optional        :: nxin, nyin
        
        ! Auxiliary
        integer(I8)                         :: nx, ny
        real(R8), allocatable               :: xgv(:), ygv(:), xg(:), &
            yg(:), vg(:), xe(:, :), ye(:, :)
        real(R8)                            :: xmin, ymin, xmax, ymax, &
            offsetx, offsety, dx, dy, dxgv, dygv

        ! Loop
        integer(I8)                         :: k

        ! Construct a 2D grid
        !====================
        ! Set mesh size
        if (present(nxin)) then 
            nx = nxin 
        else 
            nx = 200
        end if 
        if (present(nyin)) then 
            ny = nyin 
        else
            ny = 200
        end if 

        ! Allocate
        allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny), &
            vg(nx*ny))

        ! Set the gridding vectors
        call plf%ps%GetEdges(xe, ye)
        xmin = minval(xe)
        xmax = maxval(xe)
        ymin = minval(ye)
        ymax = maxval(ye)

        dx = (xmax - xmin) 
        dy = (ymax - ymin)
        
        offsetx = 0.0*dx 
        offsety = 0.0*dy 

        dxgv = (dx + 2*offsetx)/nx 
        dygv = (dy + 2*offsety)/ny 

        xgv(:) = dxgv*[(k, k = 0, nx-1)] - offsetx + xmin
        ygv(:) = dygv*[(k, k = 0, ny-1)] - offsety + ymin
        
        ! Construct
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

        ! Evaluate
        !=========
        ! Call evaluator
        call plf%Evaluate(xg, yg, 0, 0, vg)

        ! Write data
        !===========
        call Write3DCoordinateData(xg, yg, vg, savefilepath)

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(xgv, ygv, xg, yg, vg)

    end subroutine

    ! Assignment
    subroutine AssignPLF2DClass(a, b)

        class(PolygonLevelsetFunction2DUDT), allocatable, intent(inout)    :: a 
        class(PolygonLevelsetFunction2DUDT), intent(in)                    :: b 

        if (allocated(a)) then 
            deallocate(a)
        end if 

        select type (b)

        class default 

            call gdErrorHandler('Unknown type')

        type is (PolygonLevelsetFunction2DClosedApproximationUDT)

            allocate(a, source=b)
            select type (a)
            type is (PolygonLevelsetFunction2DClosedApproximationUDT)
                a = b 
            end select

        type is (PolygonLevelsetFunction2DGeneralUDT)

            allocate(a, source=b)
            select type (a)
            type is (PolygonLevelsetFunction2DGeneralUDT)
                a = b 
            end select

        type is (PolygonLevelsetFunction2DClosedExactUDT)

            allocate(a, source=b)
            select type (a)
            type is (PolygonLevelsetFunction2DClosedExactUDT)
                a = b
            end select

        end select
    
    end subroutine

    !------------------------------------------------------------------!
    !                           OPTION ROUTINES                        !
    !------------------------------------------------------------------!

    ! Option setter, general
    subroutine SetOptionsGeneral(options)

        ! Description
        !============
        ! Set options for general plf description
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DGeneralOptionsUDT)   :: options 

    end subroutine

    ! Option setter, closed exact
    subroutine SetOptionsClosedExact(options)

        ! Description
        !============
        ! Set options for general plf description
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DClosedExactOptionsUDT)   :: options 

    end subroutine

    ! Options setter, closed approximation
    subroutine SetOptionsClosedApproximation(options, resx, resy, &
        offsetx, offsety, C, M, meth)

        ! Description
        !============
        ! Set options for closed approximation plf description. One 
        ! can set the interpolant resolution (x, y direction), offset
        ! relative to extend of x and y of the polygonset, and the 
        ! continuity parameters C and M of the polygon constructor. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PLF2DClosedApproximationOptionsUDT)   :: options 
        integer(I8), intent(in)                     :: resx, resy, C, M
        real(R8), intent(in)                        :: offsetx, offsety
        character(*), intent(in)            :: meth 
        
        ! Set options
        !============
        options%resx = resx 
        options%resy = resy 
        options%offsetx = offsetx 
        options%offsety = offsety 
        options%C = C 
        options%M = M 
        options%meth = meth

        call options%optionsClosedExact%Set()

    end subroutine

    !------------------------------------------------------------------!
    !                           PLF GENERAL                            !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DGeneral(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DGeneralUDT)  :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        integer(I8)                                 :: nc

        real(R8), allocatable                       :: xe(:, :), &
            ye(:, :), nx(:), ny(:), nn(:), tx(:), ty(:), tn(:)

        ! Initialize
        !===========
        ! Add polygonset structure
        plf%ps = ps 

        ! Get edge data
        call ps%GetEdges(xe, ye)
        call ps%GetNormals(nx, ny, nn)
        call ps%GetTangents(tx, ty, tn)
        nc = size(xe, 1)
        
        ! Add to plf
        plf%x   = xe 
        plf%y   = ye
        plf%xf  = 0.5*( xe(:, 1) + xe(:, 2) )
        plf%yf  = 0.5*( ye(:, 1) + ye(:, 2) )
        plf%nxp = nx/nn ! normalize
        plf%nyp = ny/nn
        plf%txp = tx
        plf%typ = ty 
        plf%nc  = size(xe, 1)
        
    end subroutine

    ! Evaluation
    subroutine EvaluatePLF2DGeneral(plf, xq, yq, &
        derivx, derivy, vq, varin, valuesin, dplfdvarin)

        ! Description
        !============
        ! Evaluate the polygon levelset function in the query points xq, yq
        ! (value returned is in vq). The 'derivx' and 'derivy' input 
        ! arguments are used to specify the order of the derivative to 
        ! return. 

        ! Note that derivates are returned pointwise as there is no 
        ! dependence of the values on other sample points or polygon
        ! coordinates. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DGeneralUDT)  :: plf 
        real(R8), intent(in)                        :: xq(:), yq(:)
        real(R8), intent(out)                       :: vq(size(xq)) 
        integer(I8), intent(in)                     :: derivx, derivy 

        ! Optional arguments
        character(*), intent(in), optional      :: varin 
        real(R8), intent(in), optional          :: valuesin(:)
        type(MySparseUDT), optional             :: dplfdvarin
        
        character(:), allocatable               :: var 
        real(R8), allocatable                   :: values(:)
        type(MySparseUDT)                       :: dplfdvar

        ! Auxiliary
        character(:), allocatable               :: deriv, derivxc, &
            derivyc 
        integer(I8)                             :: nq, nc, nval , npsc
        integer(I8), allocatable                :: psedges(:, :), &
            edgeID(:)
        real(R8), allocatable                   :: vnx(:), vny(:), &
            vtx(:, :), vty(:, :), f(:), ind(:), fq(:), vtxq(:, :), &
            vtyq(:, :), dvnq(:), myones(:), vnxq(:), vnyq(:), dvn(:), &
            nxpq(:), nypq(:), txpq(:), typq(:), dfdxp1(:), dfdxp2(:), &
            dfdyp1(:), dfdyp2(:), xp1(:), xp2(:), yp1(:), yp2(:)
        real(R8)                                :: macheps = 0

        ! Loop
        integer(I8)                             :: i, thisloc

        ! Data
        
        ! Check input
        !============
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
        if (present(dplfdvarin)) then 
            dplfdvar = dplfdvarin 
        else
            dplfdvar = SpZeros(0, 0)
        end if 

        ! Check if derivatives are implemented
        nval = size(values, 1)
        select case(var)

        case ('no')

            ! All good

        case ('polygonsetcoordinates')

            ! Initialize basic quantities
            npsc = nval/2 ! nval should be even

            ! Initialize mapping
            call plf%ps%GetEdges(psedges) 

            ! Initialize derivative structure
            dplfdvar%ncol = nval
            dplfdvar%nrow = size(vq, 1)
            dplfdvar%nval = 4*size(vq, 1)
            call dplfdvar%Allocate()
            
        case default

            ! All bad
            call gdErrorHandler('EvaluatePLF2DGeneral: variable not implemented')

        end select

        ! Initialize
        !===========
        ! Set macheps
        macheps = 1e-13

        ! Get the number of points to evaluate on, assume proper
        ! dimensions of inputs
        nq = size(xq) 
        nc = plf%nc

        ! Allocate
        allocate(vnx(nc), vny(nc), vtx(nc, 2), vty(nc, 2), f(nc), &
            dvn(nc), ind(nc), fq(nq), vtxq(nq, 2), vtyq(nq, 2), &
            dvnq(nq), myones(nq), vnxq(nq), vnyq(nq), nxpq(nq), &
            nypq(nq), txpq(nq), typq(nq), edgeID(nq))
        
        ! Initialize
        myones = 1

        ! Associate
        associate(&
            x       => plf%x,   &
            y       => plf%y,   &
            xf      => plf%xf,  &
            yf      => plf%yf,  &
            nxp     => plf%nxp, &
            nyp     => plf%nyp, &
            txp     => plf%txp, &
            typ     => plf%typ)

        ! Compute
        !========
        ! Loop over all points
        do i = 1, nq 

            ! Normal distance
            !----------------
            ! Compute vectors between query point and edge center
            vnx = -(xf - xq(i)) 
            vny = -(yf - yq(i))

            ! Project 
            dvn = vnx * nxp + vny * nyp 

            ! Tangential distances
            !---------------------
            ! Compute vectors between query point and vertices
            vtx = x - xq(i)
            vty = y - yq(i)
            
            ! Compute 
            !--------
            ! Value
            f = abs( abs(vtx(:,1)) + abs(vtx(:, 2)) - abs(txp) ) &
                + abs( abs(vty(:,1)) + abs(vty(:, 2)) - abs(typ) ) &
                + abs(dvn)

            ! Maximum location
            thisloc = minloc(f, 1)

            ! Store variables for later
            dvnq(i)     = dvn(thisloc)
            fq(i)       = f(thisloc)
            vnxq(i)     = vnx(thisloc)
            vnyq(i)     = vny(thisloc)
            vtxq(i, :)  = vtx(thisloc, :)
            vtyq(i, :)  = vty(thisloc, :)
            txpq(i)     = txp(thisloc)
            typq(i)     = typ(thisloc)
            nxpq(i)     = nxp(thisloc)
            nypq(i)     = nyp(thisloc)
            edgeID(i)   = thisloc

        end do

        ! Hedge for machine precision
        where (fq .lt. macheps) fq = 0

        ! Return
        !=======
        ! Check what to compute exactly
        if ( (derivx > 9) .or. (derivy > 9) ) then 
            ! Throw error, these derivatives should anyway not be necessary...
            call gdErrorHandler('Derivatives are of too high order, cannot compute')
        end if
        allocate(derivxc, derivyc, deriv, source=' ')
        derivxc = ' ' 
        derivyc = ' ' 
        write(derivxc, '(I0)') derivx 
        write(derivyc, '(I0)') derivy 
        deriv = trim(derivxc) // trim(derivyc)

        select case(trim(deriv))

        case ('00')

            ! Function evaluation
            vq = fq

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Compute derivatives w.r.t. edge coordinates
                xp1 = x(edgeID, 1)
                xp2 = x(edgeID, 2)
                yp1 = y(edgeID, 1)
                yp2 = y(edgeID, 2) 
                dfdxp1 = sign(myones, abs(vtxq(:, 1)) - abs(txpq) + abs(vtxq(:, 2)))* &
                    (sign(myones, txpq) + sign(myones, vtxq(:, 1))) - &
                    sign(myones, dvnq)*(nxpq/2 + vnyq/(txpq**2 + typq**2)**(0.5) - &
                    (txpq**2*vnyq)/(txpq**2 + typq**2)**(1.5) + &
                    (txpq*typq*vnxq)/(txpq**2 + typq**2)**(1.5)) ! xp1
                dfdxp2 = - sign(myones, abs(vtxq(:, 1)) - abs(txpq) + abs(vtxq(:, 2)))* &
                    (sign(myones, txpq) - sign(myones, vtxq(:, 2))) - &
                    sign(myones, dvnq)*(nxpq/2 - vnyq/(txpq**2 + typq**2)**(0.5) + &
                    (txpq**2*vnyq)/(txpq**2 + typq**2)**(1.5) - &
                    (txpq*typq*vnxq)/(txpq**2 + typq**2)**(1.5)) ! xp2
                dfdyp1 = sign(myones, abs(vtyq(:, 1)) - abs(typq) + abs(vtyq(:, 2)))* &
                    (sign(myones, typq) + sign(myones, vtyq(:, 1))) - &
                    sign(myones, dvnq)*(nypq/2 - vnxq/(txpq**2 + typq**2)**(0.5) + &
                    (typq**2*vnxq)/(txpq**2 + typq**2)**(1.5) - &
                    (txpq*typq*vnyq)/(txpq**2 + typq**2)**(1.5)) ! yp1
                dfdyp2 = - sign(myones, abs(vtyq(:, 1)) - abs(typq) + abs(vtyq(:, 2)))* &
                    (sign(myones, typq) - sign(myones, vtyq(:, 2))) - &
                    sign(myones, dvnq)*(nypq/2 + vnxq/(txpq**2 + typq**2)**(0.5) - &
                    (typq**2*vnxq)/(txpq**2 + typq**2)**(1.5) + &
                    (txpq*typq*vnyq)/(txpq**2 + typq**2)**(1.5)) ! yp2

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [psedges(edgeID, 1), psedges(edgeID, 2), &
                    psedges(edgeID, 1)+npsc, psedges(edgeID, 2)+npsc] ! first x, then y
                dplfdvar%row = reshape(spread([(i, i = 1, nq)], 2, 4), [4*nq])
                dplfdvar%val = [dfdxp1, dfdxp2, dfdyp1, dfdyp2]

            end select

        case ('10')

            ! fx
            vq = nxpq*sign(myones, dvnq) - sign(myones, abs(vtxq(:, 1)) &
                - abs(txpq) + abs(vtxq(:, 2))) * &
                (sign(myones, vtxq(:, 1)) + sign(myones, vtxq(:, 2)))

            ! Set the gradient to 1 - can choose this value, since the 
            ! function is exactly zero here anyway. 
            where (fq .eq. 0) vq = 1

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Compute derivatives w.r.t. edge coordinates (note: dirac functions neglected)
                xp1 = x(edgeID, 1)
                xp2 = x(edgeID, 2)
                yp1 = y(edgeID, 1)
                yp2 = y(edgeID, 2) 
                dfdxp1 =  - (txpq*typq*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! xqxp1
                dfdxp2 =  (txpq*typq*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! xqxp2
                dfdyp1 = sign(myones, dvnq)/(txpq**2 + typq**2)**(0.5) - &
                    (typq**2*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! xqyp1
                dfdyp2 = (typq**2*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) - &
                    sign(myones, dvnq)/(txpq**2 + typq**2)**(0.5) ! xqyp2

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [psedges(edgeID, 1), psedges(edgeID, 2), &
                    psedges(edgeID, 1)+npsc, psedges(edgeID, 2)+npsc] ! first x, then y
                dplfdvar%row = reshape(spread([(i, i = 1, nq)], 2, 4), [4*nq])
                dplfdvar%val = [dfdxp1, dfdxp2, dfdyp1, dfdyp2]

            end select

        case ('01')

            ! fy
            vq = nypq*sign(myones, dvnq) - sign(myones, abs(vtyq(:, 1)) &
                - abs(typq) + abs(vtyq(:, 2))) * &
                (sign(myones, vtyq(:, 1)) + sign(myones, vtyq(:, 2)))

            ! Set the gradient to 1 - can choose this value, since the 
            ! function is exactly zero here anyway. 
            where (fq .eq. 0) vq = 1

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Compute derivatives w.r.t. edge coordinates (note: dirac functions neglected)
                xp1 = x(edgeID, 1)
                xp2 = x(edgeID, 2)
                yp1 = y(edgeID, 1)
                yp2 = y(edgeID, 2) 
                dfdxp1 = (txpq**2*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) &
                    - sign(myones, dvnq)/(txpq**2 + typq**2)**(0.5) ! yqxp1
                dfdxp2 = sign(myones, dvnq)/(txpq**2 + typq**2)**(0.5) - &
                    (txpq**2*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! yqxp2
                dfdyp1 = (txpq*typq*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! yqyp1
                dfdyp2 = - (txpq*typq*sign(myones, dvnq))/(txpq**2 + typq**2)**(1.5) ! yqyp2

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [psedges(edgeID, 1), psedges(edgeID, 2), &
                    psedges(edgeID, 1)+npsc, psedges(edgeID, 2)+npsc] ! first x, then y
                dplfdvar%row = reshape(spread([(i, i = 1, nq)], 2, 4), [4*nq])
                dplfdvar%val = [dfdxp1, dfdxp2, dfdyp1, dfdyp2]

            end select

        case ('20')

            ! fxx - simply zero
            vq(:) = 0
            dplfdvar = SpZeros(nq, nval)

        case ('11')

            ! fxy - simply zero
            vq(:) = 0
            dplfdvar = SpZeros(nq, nval)

        case ('02')

            ! fyy - simply zero
            vq(:) = 0
            dplfdvar = SpZeros(nq, nval)

        case default

            ! Put to zero - should be higher order derivatve
            vq(:) = 0

        end select

        ! Housekeeping
        !=============
        ! Optional output arguments
        if (present(dplfdvarin)) then 
            dplfdvarin = dplfdvar
        end if 

        ! End associate
        end associate

        ! Deallocate
        deallocate(vnx, vny, vtx, vty, f, ind, fq, vtxq, vtyq, dvnq)

    end subroutine

    !------------------------------------------------------------------!
    !                         PLF CLOSED EXACT                         !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DClosedExact(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedExactUDT)      :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        integer(I8)                                 :: np, flag, ne, nl
        integer(I8), allocatable                    :: tempvertlabel(:, :), &
            tempedgelabel(:, :), vertlabelp(:, :), edgelabelp(:, :), &
            vp(:), vp1(:), vp2(:)

        real(R8), allocatable                       :: xe(:, :), &
            ye(:, :), nx(:), ny(:), nn(:), tx(:), ty(:), tn(:), xp(:), &
            yp(:), tempx(:), tempy(:), tempnx(:, :), tempny(:, :), &
            nxpv(:, :), nypv(:, :), theta0(:), crossprod(:), tnxp(:), &
            tnyp(:), tnn(:), xf(:), yf(:), xp1(:), yp1(:), xp2(:), &
            yp2(:), txp1(:), typ1(:), txp2(:), typ2(:)

        ! Loop
        integer(I8)                                 :: i, ce, ce2, k 
        ! Initialize
        !===========
        ! Check polygonset 
        if (ps%np > 0) then ! allow for empty polygon sets - will evaluate to zero
            call ps%OrientNestedClosedPolygons(flag)
            if (flag > 1) then 
                ! Polygon set does not comply, throw error
                call gdErrorHandler('InitializePLF2DClosedExact: polygon set should be closed and non-intersecting')
            end if 
        else
            ! Allocate all to zero sized arrays
            allocate(plf%xp(0), plf%yp(0), plf%xf(0), plf%yf(0), &
                plf%nxp(0), plf%nyp(0), plf%tnp(0), plf%nxpv(0, 2), &
                plf%nypv(0, 2), plf%crossprod(0), plf%thetav(0), plf%xp1(0), &
                plf%xp2(0), plf%yp1(0), plf%yp2(0))
            plf%ps = ps 
            return 
        end if 
        
        ! Add polygonset structure
        plf%ps = ps 

        ! Associate
        associate(&
            pol     => ps%polygons &
        )

        ! Compute additional data
        !========================
        ! Compute number of labels
        nl = size(pol(1)%labels, 2)

        ! Get edge data
        call ps%GetEdges(xe, ye) ! may still change
        call ps%GetNormals(nx, ny, nn) ! may still change
        call ps%GetTangents(tx, ty, tn)
        !call ps%GetPoints(xp, yp) ! may still change

        ! Allocate
        np = (size(xe, 1))+ps%np
        allocate(nxpv(np, 2), nypv(np, 2), xp(np), yp(np), &
            xf(size(xe, 1)), yf(size(xe, 1)), xp1(size(xe, 1)),  &
            yp1(size(xe, 1)), xp2(size(xe, 1)), yp2(size(xe, 1)), &
            vertlabelp(np, nl), edgelabelp(size(xe, 1), nl), &
            vp(np), vp1(size(xe, 1)), vp2(size(xe, 1)))

        ! Compute vertex regions
        ce = 0
        ce2 = 0
        do i = 1, ps%np 
            ! Normally, since each polygon is closed and non-intersecting,
            ! the first and last vertex will be the same and no 
            ! duplicate vertices will exist in vert. 
            
            ! Unpack
            ne = pol(i)%ne 

            ! Coordinates
            tempx = pol(i)%x(pol(i)%vert)
            tempy = pol(i)%y(pol(i)%vert)
            tempvertlabel = pol(i)%labels(pol(i)%vert, :)

            txp1 = tempx(1:ne)
            txp2 = tempx(2:ne+1)
            typ1 = tempy(1:ne)
            typ2 = tempy(2:ne+1)
            tempedgelabel = tempvertlabel(1:ne, :) ! somewhat extrapolation

            ! Normals
            allocate(tempnx(ne+1, 2), tempny(ne+1, 2))

            tnxp = -(typ2 - typ1)
            tnyp = (txp2 - txp1)
            tnn = sqrt(tnxp**2 + tnyp**2)
            tnxp = tnxp/tnn 
            tnyp = tnyp/tnn

            tempnx(:, 1) = [tnxp(ne), tnxp]
            tempnx(:, 2) = [tnxp, tnxp(1)]
            tempny(:, 1) = [tnyp(ne), tnyp]
            tempny(:, 2) = [tnyp, tnyp(1)]



            !tempnx(1:ne, 2) = pol(i)%nx/pol(i)%nn 
            !tempnx(1:ne, 1) = [pol(i)%nx(ne)/pol(i)%nn(ne), &
            !    pol(i)%nx(1:ne-1)/pol(i)%nn(1:ne-1)]
            !tempnx(ne+1, 1:2) = [pol(i)%nx(1)/pol(i)%nn(1), pol(i)%nx(ne)/pol(i)%nn(ne)]

            !tempny(1:ne, 2) = pol(i)%ny/pol(i)%nn 
            !tempny(1:ne, 1) = [pol(i)%ny(ne)/pol(i)%nn(ne), &
            !    pol(i)%ny(1:ne-1)/pol(i)%nn(1:ne-1)]
            !tempny(ne+1, 1:2) = [pol(i)%ny(1)/pol(i)%nn(1), pol(i)%ny(ne)/pol(i)%nn(ne)]


            ! Add
            xp(ce+1:ce+ne+1) = tempx 
            yp(ce+1:ce+ne+1) = tempy 
            vp(ce+1:ce+ne+1) = [(k, k = ce+1, ce+ne+1)]
            vertlabelp(ce+1:ce+ne+1, :) = tempvertlabel
            edgelabelp(ce2+1:ce2+ne, :) = tempedgelabel
            nxpv(ce+1:ce+ne+1, :) = tempnx ! normals should be oriented outwards already
            nypv(ce+1:ce+ne+1, :) = tempny 
            nx(ce2+1:ce2+ne) = tnxp !pol(i)%nx/pol(i)%nn 
            ny(ce2+1:ce2+ne) = tnyp! pol(i)%ny/pol(i)%nn 
            tn(ce2+1:ce2+ne) = pol(i)%nn
            xf(ce2+1:ce2+ne) = 0.5*(tempx(1:ne) + tempx(2:ne+1)) 
            yf(ce2+1:ce2+ne) = 0.5*(tempy(1:ne) + tempy(2:ne+1)) 
            xp1(ce2+1:ce2+ne) = txp1
            xp2(ce2+1:ce2+ne) = txp2
            yp1(ce2+1:ce2+ne) = typ1
            yp2(ce2+1:ce2+ne) = typ2
            vp1(ce2+1:ce2+ne) = vp(ce+1:ce+ne)
            vp2(ce2+1:ce2+ne) = vp(ce+2:ce+ne+1)


            ! Update counter
            ce = ce + ne + 1
            ce2 = ce2 + ne

            ! Housekeeping
            deallocate(tempnx, tempny)

            
            
        end do 

        ! Compute cross product
        crossprod = nxpv(:, 2)*nypv(:, 1) - nxpv(:, 1)*nypv(:, 2)

        ! Compute angles
        theta0 = atan2(nxpv(:, 2)*nypv(:, 1) - nxpv(:, 1)*nypv(:, 2), &
            nxpv(:, 1)*nxpv(:, 2) + nypv(:, 1)*nypv(:, 2))
        where (theta0 < 0) theta0 = theta0 + 2*pi_R8         
        
        ! Add to plf
        plf%xp      = xp 
        plf%yp      = yp
        plf%vp      = vp 
        plf%xp1     = xp1 
        plf%yp1     = yp1
        plf%xp2     = xp2 
        plf%yp2     = yp2
        plf%xf      = xf
        plf%yf      = yf
        plf%nxp     = nx 
        plf%nyp     = ny
        plf%tnp     = tn 
        plf%nxpv    = nxpv 
        plf%nypv    = nypv
        plf%thetav  = theta0 
        plf%crossprod   = crossprod
        plf%edgelabel = edgelabelp 
        plf%vertlabel = vertlabelp
        plf%vp1     = vp1
        plf%vp2     = vp2 

        ! Housekeeping
        !=============
        end associate

    end subroutine 

    ! Evaluation
    subroutine EvaluatePLF2DClosedExact(plf, xq, yq, derivx, derivy, vq, &
        varin, valuesin, dplfdvarin)

        ! Description
        !============
        ! Evaluate the polygon levelset function in the query points xq, yq
        ! (value returned is in vq). The 'derivx' and 'derivy' input 
        ! arguments are used to specify the order of the derivative to 
        ! return. 

        ! Note that derivates are returned pointwise as there is no 
        ! dependence of the values on other sample points or polygon
        ! coordinates.
        
        use iso_fortran_env
        use, intrinsic :: ieee_arithmetic

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedExactUDT)  :: plf 
        real(R8), intent(in)                        :: xq(:), yq(:)
        real(R8), intent(out)                       :: vq(size(xq)) 
        integer(I8), intent(in)                     :: derivx, derivy 

        ! Optional arguments
        character(*), intent(in), optional      :: varin 
        real(R8), intent(in), optional          :: valuesin(:)
        type(MySparseUDT), optional             :: dplfdvarin

        character(:), allocatable               :: var 
        real(R8), allocatable                   :: values(:)
        type(MySparseUDT)                       :: dplfdvar

        ! Auxiliary
        character(:), allocatable               :: deriv, derivxc, &
            derivyc

        integer(I8)                             :: nq, np, npe, &
            indmine, indminv, indmin, nval, ne, npsc
        integer(I8), allocatable                :: eind(:), vind(:), &
            minind(:), psvert(:), psedges(:, :)

        real(R8)                                :: inf, signe, signv, &
            fv, fe, xqr, yqr, totsign(1:2)
        real(R8), parameter                     :: myone = 1
        real(R8), allocatable                   :: myones(:), dvn(:), &
            tdistvert(:), tcrossprod(:), vx(:), vy(:), tvn(:), &
            dx(:), dy(:), theta(:), distedge(:), distvert(:), &
            val(:), tvx(:), tvy(:), ttxp(:), ttyp(:), tnxp(:), &
            tnyp(:), dfdxp1(:), dfdxp2(:), dfdyp1(:), dfdyp2(:), &
            dfdxp(:), dfdyp(:)
        

        logical, allocatable                    :: isinvert(:), &
            onedge(:), te(:), tv(:)

        ! Loop
        integer(I8)                             :: i, iq, ccv, cce, cnv

        ! Data
        inf = ieee_value(inf, ieee_positive_inf)

        
        ! Check input
        !============
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
        if (present(dplfdvarin)) then 
            dplfdvar = dplfdvarin 
        end if 

        if (plf%ps%np == 0) then 
            ! Empty polygonset - return zero values
            val = 0 
            return 
        end if 

        ! Initialize
        !===========
        ! Associate
        associate(&
            xp      => plf%xp,      &
            yp      => plf%yp,      &
            xp1     => plf%xp1,     &
            xp2     => plf%xp2,     &
            yp1     => plf%yp1,     &
            yp2     => plf%yp2,     &
            nxp     => plf%nxp,     &
            nyp     => plf%nyp,     &
            nxpv    => plf%nxpv,    &
            nypv    => plf%nypv,    &
            tnp     => plf%tnp,     &
            theta0  => plf%thetav,  &
            crossprod => plf%crossprod, &
            xf      => plf%xf,      &
            yf      => plf%yf       &
            )

        ! Get sizes
        npe = size(xp1, 1)
        np = size(xp, 1)
        nq = size(xq, 1)

        ! Check if derivatives are implemented
        nval = size(values, 1)
        select case(var)

        case ('no')

            ! All good

        case ('polygonsetcoordinates')

            ! Initialize basic quantities
            npsc = nval/2 ! nval should be even

            ! Initialize mapping
            allocate(psvert(np), psedges(npe, 2))
            ccv = 0
            cce = 0
            cnv = 0
            do i = 1, plf%ps%np 
                ! Add
                ne = plf%ps%polygons(i)%ne 
                psvert(ccv+1:ccv+ne+1) = plf%ps%polygons(i)%vert + cnv
                psedges(cce+1:cce+ne, 1) = plf%ps%polygons(i)%vert(1:ne) + cnv 
                psedges(cce+1:cce+ne, 2) = plf%ps%polygons(i)%vert(2:ne+1) + cnv

                ! Update counters
                ccv = ccv + ne + 1
                cce = cce + ne
                cnv = cnv + plf%ps%polygons(i)%nv
            end do

            ! Initialize derivative structure (note: number values 
            ! depends, initialized later)
            dplfdvar%ncol = nval
            dplfdvar%nrow = size(vq, 1)
            
        case default

            ! All bad
            call gdErrorHandler('EvaluatePLF2DGeneral: variable not implemented')

        end select
        
        ! Allocate and initialize
        allocate(eind(nq), vind(nq), minind(nq), tdistvert(nq), &
            tcrossprod(nq), val(nq))

        ! Check what to compute exactly
        if ( (derivx > 9) .or. (derivy > 9) ) then 
            ! Throw error, these derivatives should anyway not be necessary...
            call gdErrorHandler('Derivatives are of too high order, cannot compute')
        end if
        allocate(derivxc, derivyc, deriv, source=' ')
        derivxc = ' ' 
        derivyc = ' ' 
        write(derivxc, '(I0)') derivx 
        write(derivyc, '(I0)') derivy 
        deriv = trim(derivxc) // trim(derivyc)

        ! Loop 
        !=====
        allocate(myones(np))
        myones = 1
        !$omp parallel do default(none) shared(xq, yq, plf, &
        !$omp val, minind, vind, eind, tdistvert, tcrossprod, nq, myones, inf) &
        !$omp private(xqr, yqr, vx, vy, dvn, tvn, dx, dy, theta, isinvert, &
        !$omp onedge, distedge, distvert, indmine, fe, signe, indminv, &
        !$omp fv, signv, totsign, indmin) if(.not. omp_in_parallel()) &
        !$omp schedule(static)
        do iq = 1, nq
            
            ! Unpack
            xqr = xq(iq)
            yqr = yq(iq)

            ! Construct vector between edge centers and points
            vx = (xqr - xf)
            vy = (yqr - yf)
            dvn = (vx*nxp + vy*nyp)
            tvn = vx*nyp - vy*nxp ! distance in tangential direction
            
            ! Check in which regions the query points lie
            ! Vertex regions
            dx = xqr - xp
            dy = yqr - yp
            theta = atan2(dx*nypv(:, 1) - dy*nxpv(:, 1), dx*nxpv(:, 1) + dy*nypv(:, 1))
            where (theta < 0) theta = theta + 2*pi_R8
            isinvert = theta <= theta0

            ! Hedge for vertices lying exactly on polygonset vertex
            where ((dx == 0) .and. (dy == 0)) isinvert = .true. 
            
            ! Edge regions
            onedge = (abs(tvn) <= 0.5*tnp)
            
            ! Compute distances
            distedge = dvn
            distvert = sign(myones, crossprod)*sqrt((xp - xqr)**2 + (yp - yqr)**2)
            where (.not. onedge) distedge = inf 
            where (.not. isinvert) distvert = inf 
            
            ! Compute minimal distance
            indmine = minloc(abs(distedge), 1)
            fe = minval(abs(distedge))
            signe = sign(myone, distedge(indmine))
            eind(iq) = indmine

            indminv = minloc(abs(distvert), 1)
            fv = minval(abs(distvert))
            tdistvert(iq) = abs(distvert(indminv))
            tcrossprod(iq) = crossprod(indminv)

            signv = sign(myone, distvert(indminv))
            vind(iq) = indminv
            totsign = [signe, signv]
            val(iq) = minval([fe, fv])
            indmin = minloc([fe, fv], 1)
            minind(iq) = indmin 
            val(iq) = val(iq)*totsign(indmin)
            
        end do 
        !$omp end parallel do

        ! Construct function handle depending on type
        select case (trim(deriv))
        
        case ('00')
            
            ! Function value
            vq = val
            
            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Precompute values
                tvx = (xq - xf(eind))
                tvy = (yq - yf(eind))
                ttxp = (xp2(eind) - xp1(eind))
                ttyp = (yp2(eind) - yp1(eind))
                tnxp = nxp(eind)
                tnyp = nyp(eind)
                if (allocated(myones)) then 
                    deallocate(myones)
                end if 
                allocate(myones(nq))
                myones = 1
            
                te = minind == 1
                tv = minind == 2

                ! Compute derivatives w.r.t. edge coordinates
                allocate(dfdxp1(nq), dfdxp2(nq), dfdyp1(nq), &
                    dfdyp2(nq), dfdxp(nq), dfdyp(nq))
                where (te)
                    dfdxp1 = (ttxp**2*tvy)/(ttxp**2 + ttyp**2)**(1.5) - &
                        tvy/(ttxp**2 + ttyp**2)**(0.5) - tnxp*0.5 - &
                        (ttxp*ttyp*tvx)/(ttxp**2 + ttyp**2)**(1.5) ! xp1
                    dfdxp2 = tvy/(ttxp**2 + ttyp**2)**(0.5) - tnxp*0.5 - &
                        (ttxp**2*tvy)/(ttxp**2 + ttyp**2)**(1.5) + &
                        (ttxp*ttyp*tvx)/(ttxp**2 + ttyp**2)**(1.5) ! xp2
                    dfdyp1 = tvx/(ttxp**2 + ttyp**2)**(0.5) - tnyp*0.5 - &
                        (ttyp**2*tvx)/(ttxp**2 + ttyp**2)**(1.5) + &
                        (ttxp*ttyp*tvy)/(ttxp**2 + ttyp**2)**(1.5) ! yp1
                    dfdyp2 = (ttyp**2*tvx)/(ttxp**2 + ttyp**2)**(1.5) - &
                        tvx/(ttxp**2 + ttyp**2)**(0.5) - tnyp*0.5 - &
                        (ttxp*ttyp*tvy)/(ttxp**2 + ttyp**2)**(1.5) ! yp2
                elsewhere
                    dfdxp1 = 0
                    dfdyp1 = 0
                    dfdxp2 = 0
                    dfdyp2 = 0
                end where

                where (tv)
                    dfdxp = sign(myones, tcrossprod)*(xp(vind) - xq)/(abs(tdistvert))
                    dfdyp = sign(myones, tcrossprod)*(yp(vind) - yq)/(abs(tdistvert))
                elsewhere 
                    dfdxp = 0
                    dfdyp = 0
                end where

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [pack(psedges(eind, 1), te), pack(psedges(eind, 2), te), &
                    pack(psedges(eind, 1), te)+npsc, pack(psedges(eind, 2), te)+npsc, &
                    pack(psvert(vind), tv), pack(psvert(vind), tv)+npsc]
                dplfdvar%row = [reshape(spread(pack([(i, i = 1, nq)], te), 2, 4), [4*count(te)]), &
                    reshape(spread(pack([(i, i = 1, nq)], tv), 2, 4), [2*count(tv)])]
                dplfdvar%val = [pack(dfdxp1, te), pack(dfdxp2, te), &
                    pack(dfdyp1, te), pack(dfdyp2, te), &
                    pack(dfdxp, tv), pack(dfdyp, tv)]
                dplfdvar%nval = size(dplfdvar%val, 1)

            end select

        case ('10')
            
            ! First order derivative w.r.t. xq
            
            ! Compute
            te = minind == 1
            tv = minind == 2
            
            do iq = 1, nq
                if (te(iq)) then 
                    vq(iq) = nxp(eind(iq))
                else
                    if (tdistvert(iq) == 0.0) then 
                        vq(iq) = 1
                    else
                        vq(iq) = -sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)*(xp(vind(iq)) - xq(iq))
                    end if 
                end if
            end do 

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Precompute values
                tvx = (xq - xf(eind))
                tvy = (yq - yf(eind))
                ttxp = (xp2(eind) - xp1(eind))
                ttyp = (yp2(eind) - yp1(eind))
                tnxp = nxp(eind)
                tnyp = nyp(eind)
                if (allocated(myones)) then 
                    deallocate(myones)
                end if 
                allocate(myones(nq))
                myones = 1

                te = minind == 1
                tv = minind == 2

                ! Compute derivatives w.r.t. edge coordinates
                allocate(dfdxp1(nq), dfdxp2(nq), dfdyp1(nq), &
                    dfdyp2(nq), dfdxp(nq), dfdyp(nq))
                where (te)
                    dfdxp1 = -(ttxp*ttyp)/(ttxp**2 + ttyp**2)**(1.5) ! xqxp1
                    dfdxp2 = (ttxp*ttyp)/(ttxp**2 + ttyp**2)**(1.5) ! xqxp2
                    dfdyp1 = 1.0/(ttxp**2 + ttyp**2)**(0.5) - ttyp**2/(ttxp**2 + ttyp**2)**(1.5) ! xqyp1
                    dfdyp2 = ttyp**2/(ttxp**2 + ttyp**2)**(1.5) - 1.0/(ttxp**2 + ttyp**2)**(0.5) ! xqyp2
                elsewhere
                    dfdxp1 = 0
                    dfdyp1 = 0
                    dfdxp2 = 0
                    dfdyp2 = 0
                end where

                where (tv .and. (tdistvert /= 0.0))
                    dfdxp = (sign(myones, tcrossprod)*(xp(vind) - xq)**2)/(2*tdistvert**3) - sign(myones, tcrossprod)/tdistvert ! xqxp
                    dfdyp = (sign(myones, tcrossprod)*(xp(vind) - xq)*(yp(vind) - yq))/(tdistvert**3) ! xqyp
                elsewhere 
                    dfdxp = 0
                    dfdyp = 0
                end where

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [pack(psedges(eind, 1), te), pack(psedges(eind, 2), te), &
                    pack(psedges(eind, 1), te)+npsc, pack(psedges(eind, 2), te)+npsc, &
                    pack(psvert(vind), tv), pack(psvert(vind), tv)+npsc]
                    dplfdvar%row = [reshape(spread(pack([(i, i = 1, nq)], te), 2, 4), [4*count(te)]), &
                    reshape(spread(pack([(i, i = 1, nq)], tv), 2, 4), [2*count(tv)])]
                dplfdvar%val = [pack(dfdxp1, te), pack(dfdxp2, te), &
                    pack(dfdyp1, te), pack(dfdyp2, te), &
                    pack(dfdxp, tv), pack(dfdyp, tv)]
                dplfdvar%nval = size(dplfdvar%val, 1)

            end select
        
        case ('01')
            
            ! First order derivative w.r.t. yq
            
            ! Compute
            te = minind == 1
            tv = minind == 2

            do iq = 1, nq
                if (te(iq)) then 
                    vq(iq) = nyp(eind(iq))
                else
                    if (tdistvert(iq) == 0.0) then 
                        vq(iq) = 1
                    else
                        vq(iq) = -sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)*(yp(vind(iq)) - yq(iq))
                    end if
                end if
            end do 

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case ('polygonsetcoordinates')

                ! Derivatives w.r.t. polygon set coordinates

                ! Precompute values
                tvx = (xq - xf(eind))
                tvy = (yq - yf(eind))
                ttxp = (xp2(eind) - xp1(eind))
                ttyp = (yp2(eind) - yp1(eind))
                tnxp = nxp(eind)
                tnyp = nyp(eind)
                if (allocated(myones)) then 
                    deallocate(myones)
                end if 
                allocate(myones(nq))
                myones = 1

                te = minind == 1
                tv = minind == 2

                ! Compute derivatives w.r.t. edge coordinates
                allocate(dfdxp1(nq), dfdxp2(nq), dfdyp1(nq), &
                    dfdyp2(nq), dfdxp(nq), dfdyp(nq))
                where (te)
                    dfdxp1 = ttxp**2/(ttxp**2 + ttyp**2)**(1.5) - 1.0/(ttxp**2 + ttyp**2)**(0.5) ! yqxp1
                    dfdxp2 = 1.0/(ttxp**2 + ttyp**2)**(0.5) - ttxp**2/(ttxp**2 + ttyp**2)**(1.5) ! yqxp2
                    dfdyp1 = (ttxp*ttyp)/(ttxp**2 + ttyp**2)**(1.5) ! yqyp1
                    dfdyp2 = -(ttxp*ttyp)/(ttxp**2 + ttyp**2)**(1.5) ! yqyp2
                elsewhere
                    dfdxp1 = 0
                    dfdyp1 = 0
                    dfdxp2 = 0
                    dfdyp2 = 0
                end where

                where (tv .and. (tdistvert /= 0.0))
                    dfdxp = (sign(myones, tcrossprod)*(xp(vind) - xq)*(yp(vind) - yq))/(tdistvert**3) ! yqxp
                    dfdyp = (sign(myones, tcrossprod)*(yp(vind) - yq)**2)/(2*tdistvert**3) - sign(myones, tcrossprod)/tdistvert ! yqyp
                elsewhere 
                    dfdxp = 0
                    dfdyp = 0
                end where

                ! Compute derivatives w.r.t. polygon set coordinates
                dplfdvar%col = [pack(psedges(eind, 1), te), pack(psedges(eind, 2), te), &
                    pack(psedges(eind, 1), te)+npsc, pack(psedges(eind, 2), te)+npsc, &
                    pack(psvert(vind), tv), pack(psvert(vind), tv)+npsc]
                dplfdvar%row = [reshape(spread(pack([(i, i = 1, nq)], te), 2, 4), [4*count(te)]), &
                    reshape(spread(pack([(i, i = 1, nq)], tv), 2, 4), [2*count(tv)])]
                dplfdvar%val = [pack(dfdxp1, te), pack(dfdxp2, te), &
                    pack(dfdyp1, te), pack(dfdyp2, te), &
                    pack(dfdxp, tv), pack(dfdyp, tv)]
                dplfdvar%nval = size(dplfdvar%val, 1)

            end select

        case ('20')

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case default

                ! Not implemented
                call gdErrorHandler('EvaluatePLF2DClosedExact: derivatives not implemented')

            end select
            
            ! Second order derivative w.r.t. xq, xq
            vq = 0
            te = minind == 1
            tv = minind == 2
            do iq = 1, nq
                if (tv(iq) .and. (tdistvert(iq) /= 0.0)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(xp(vind(iq)) - xq(iq))**2 + &
                    sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)
                end if
            end do 

        case ('11')
            
            ! Second order derivative w.r.t. xq, yq
            vq = 0
            te = minind == 1;
            tv = minind == 2;
            do iq = 1, nq
                if (tv(iq) .and. (tdistvert(iq) /= 0.0)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(xp(vind(iq)) - xq(iq))*(yp(vind(iq)) - xq(iq))
                end if
            end do 

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case default

                ! Not implemented
                call gdErrorHandler('EvaluatePLF2DClosedExact: derivatives not implemented')

            end select
            
        case ('02')
            
            vq = 0
            te = minind == 1;
            tv = minind == 2;

            do iq = 1, nq
                if (tv(iq) .and. (tdistvert(iq) /= 0.0)) then 
                    vq(iq) = -sign(myone, tdistvert(iq))*&
                    sign(myone, tcrossprod(iq))/tdistvert(iq)**3&
                    *(yp(vind(iq)) - yq(iq))**2 + &
                    sign(myone, tdistvert(iq))*sign(myone, tcrossprod(iq))/tdistvert(iq)
                end if
            end do 

            ! Check derivatives
            select case(var)

            case ('no')
                
                ! No derivatives
                dplfdvar = SpZeros(nq, size(values, 1))

            case default

                ! Not implemented
                call gdErrorHandler('EvaluatePLF2DClosedExact: derivatives not implemented')

            end select
            
        case default

            ! Just for now
            vq = 0

        end select

        ! Housekeeping
        !=============
        end associate
        
        ! Optional arguments
        if (present(dplfdvarin)) then 
            dplfdvarin = dplfdvar 
        end if

    end subroutine

    ! Label evaluation
    subroutine EvaluatePLF2DClosedExactLabel(plf, xq, yq, vq, edgeIDopt, vertIDopt)

        ! Description
        !============
        ! Evaluate the polygon levelset labels in the query points xq, yq
        ! (value returned is in vq). The label is of the same size as 
        ! the source polygon labels and evaluation is done by simply 
        ! taking the corresponding polygon edge's label (or vertex) 
        ! as originally given during the construction phase. No 
        ! derivatives etc are returned. 

        ! Note: optionally, the edge/vertex ID of the PLF which was used
        ! to evaluate the label is returned in edgeIDopt/vertIDopt 
        ! (zero if not on edge/vertex)

        use iso_fortran_env
        use, intrinsic :: ieee_arithmetic

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedExactUDT)  :: plf 
        real(R8), intent(in)                        :: xq(:), yq(:)
        integer(I8), allocatable, intent(out)       :: vq(:, :) 
        integer(I8), allocatable, optional, intent(out)     :: edgeIDopt(:), &
            vertIDopt(:)

        ! Auxiliary
        integer(I8)                             :: nq, np, npe, &
            indmine, indminv, indmin, nl
        integer(I8), allocatable                :: eind(:), vind(:), &
            minind(:), edgeID(:), vertID(:)

        real(R8)                                :: inf, signe, signv, &
            fv, fe, xqr, yqr, totsign(1:2)
        real(R8), parameter                     :: myone = 1
        real(R8), allocatable                   :: myones(:), dvn(:), &
            tdistvert(:), tcrossprod(:), vx(:), vy(:), tvn(:), &
            dx(:), dy(:), theta(:), distedge(:), distvert(:), &
            val(:)

        logical, allocatable                    :: isinvert(:), &
            onedge(:)

        ! Loop
        integer(I8)                             :: iq

        ! Data
        inf = ieee_value(inf, ieee_positive_inf)        

        ! Initialize
        !===========
        ! Check for trivial case
        if (plf%ps%np == 0) then 
            ! Empty polygonset - return zero values
            vq = 0 
            if (present(edgeIDopt)) edgeIDopt = 0
            if (present(vertIDopt)) vertIDopt = 0
            return 
        end if 

        ! Associate
        associate(&
            vl      => plf%vertlabel,   &
            el      => plf%edgelabel,   &
            xp      => plf%xp,      &
            yp      => plf%yp,      &
            xp1     => plf%xp1,     &
            xp2     => plf%xp2,     &
            yp1     => plf%yp1,     &
            yp2     => plf%yp2,     &
            nxp     => plf%nxp,     &
            nyp     => plf%nyp,     &
            nxpv    => plf%nxpv,    &
            nypv    => plf%nypv,    &
            tnp     => plf%tnp,     &
            theta0  => plf%thetav,  &
            crossprod => plf%crossprod, &
            xf      => plf%xf,      &
            yf      => plf%yf       &
            )

        ! Get sizes
        npe = size(xp1, 1)
        np = size(xp, 1)
        nq = size(xq, 1)
        nl = size(vl, 2)

        ! Initialize
        allocate(vq(nq, nl))

        ! Check for the trivial case
        if (nl == 0) then 
            return 
        end if 
        
        ! Allocate and initialize
        allocate(eind(nq), vind(nq), minind(nq), tdistvert(nq), &
            tcrossprod(nq), val(nq), edgeID(nq), vertID(nq))
        edgeID = 0
        vertID = 0

        ! Loop 
        !=====
        allocate(myones(np))
        myones = 1
        !$omp parallel do default(none) shared(xq, yq, plf, vertID, edgeID, &
        !$omp val, minind, vind, eind, tdistvert, tcrossprod, nq, myones, inf, vq) &
        !$omp private(xqr, yqr, vx, vy, dvn, tvn, dx, dy, theta, isinvert, &
        !$omp onedge, distedge, distvert, indmine, fe, signe, indminv, &
        !$omp fv, signv, totsign, indmin) if (.not. omp_in_parallel()) & 
        !$omp schedule(static)
        do iq = 1, nq
            
            ! Unpack
            xqr = xq(iq)
            yqr = yq(iq)

            ! Construct vector between edge centers and points
            vx = (xqr - xf)
            vy = (yqr - yf)
            dvn = (vx*nxp + vy*nyp)
            tvn = vx*nyp - vy*nxp ! distance in tangential direction
            
            ! Check in which regions the query points lie
            ! Vertex regions
            dx = xqr - xp
            dy = yqr - yp
            theta = atan2(dx*nypv(:, 1) - dy*nxpv(:, 1), dx*nxpv(:, 1) + dy*nypv(:, 1))
            where (theta < 0) theta = theta + 2*pi_R8
            isinvert = theta <= theta0

            ! Hedge for vertices lying exactly on polygonset vertex
            where ((dx == 0) .and. (dy == 0)) isinvert = .true. 
            
            ! Edge regions
            onedge = (abs(tvn) <= 0.5*tnp)
            
            ! Compute distances
            distedge = dvn
            distvert = sign(myones, crossprod)*sqrt((xp - xqr)**2 + (yp - yqr)**2)
            where (.not. onedge) distedge = inf 
            where (.not. isinvert) distvert = inf 
            
            ! Compute minimal distance
            indmine = minloc(abs(distedge), 1)
            fe = minval(abs(distedge))
            signe = sign(myone, distedge(indmine))
            eind(iq) = indmine

            indminv = minloc(abs(distvert), 1)
            fv = minval(abs(distvert))
            tdistvert(iq) = abs(distvert(indminv))
            tcrossprod(iq) = crossprod(indminv)

            signv = sign(myone, distvert(indminv))
            vind(iq) = indminv
            totsign = [signe, signv]
            val(iq) = minval([fe, fv])
            indmin = minloc([fe, fv], 1)
            if (fv == 0) then ! Force vertex value if vertex distance is precisely zero
                indmin = 2
            end if
            minind(iq) = indmin 
            val(iq) = val(iq)*totsign(indmin)
            if (indmin == 1) then 
                vq(iq, :) = el(indmine, :)
                edgeID(iq) = indmine 
            else
                vq(iq, :) = vl(indminv, :)
                vertID(iq) = indminv 
            end if 
            
        end do 
        !$omp end parallel do


        ! Housekeeping
        !=============
        if (present(edgeIDopt)) then 
            edgeIDopt = edgeID 
        end if 
        if (present(vertIDopt)) then 
            vertIDopt = vertID 
        end if 
        end associate

    end subroutine

    ! Visualization
    subroutine VisualizePolygonLevelsetFunctionLabelClosedExact(plf, savefilepath, &
        nxin, nyin, labelindin)

    ! Description
    !============
    ! This routine provides a crude way of visualizing the 
    ! levelset function by making a 2D patchplot on a structured 
    ! grid. The extent is equal to the coordinate boundaries 
    ! of the polygon set vertices, and the resolution in both 
    ! directions is fixed and constant. 

    ! Declare variables
    !==================
    ! Arguments
    class(PolygonLevelsetFunction2DClosedExactUDT)      :: plf 
    character(*), intent(in)                 :: savefilepath
    integer(I8), intent(in), optional        :: nxin, nyin, labelindin
    
    ! Auxiliary
    integer(I8)                         :: nx, ny, labelind
    integer(I8), allocatable            :: vg(:, :)
    real(R8), allocatable               :: xgv(:), ygv(:), xg(:), &
        yg(:), xe(:, :), ye(:, :), vgr(:, :)
    real(R8)                            :: xmin, ymin, xmax, ymax, &
        offsetx, offsety, dx, dy, dxgv, dygv

    ! Loop
    integer(I8)                         :: k

    ! Construct a 2D grid
    !====================
    ! Check
    if (plf%ps%np == 0) then 
        print *,'Empty polygon levelset, cannot evaluate. Returning...'
        return 
    end if 

    ! Set mesh size
    if (present(nxin)) then 
        nx = nxin 
    else 
        nx = 200
    end if 
    if (present(nyin)) then 
        ny = nyin 
    else
        ny = 200
    end if 
    if (present(labelindin)) then 
        labelind = labelindin
    else
        labelind = 1
    end if 

    ! Allocate
    allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny))

    ! Set the gridding vectors
    call plf%ps%GetEdges(xe, ye)
    xmin = minval(xe)
    xmax = maxval(xe)
    ymin = minval(ye)
    ymax = maxval(ye)

    dx = (xmax - xmin) 
    dy = (ymax - ymin)
    
    offsetx = 0.0*dx 
    offsety = 0.0*dy 

    dxgv = (dx + 2*offsetx)/nx 
    dygv = (dy + 2*offsety)/ny 

    xgv(:) = dxgv*[(k, k = 0, nx-1)] - offsetx + xmin
    ygv(:) = dygv*[(k, k = 0, ny-1)] - offsety + ymin
    
    ! Construct
    call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

    ! Evaluate
    !=========
    ! Call evaluator
    call plf%EvaluateLabel(xg, yg, vg)

    ! Write data
    !===========
    vgr = vg*1.0_R8
    call Write3DCoordinateData(xg, yg, vgr(:, labelind), savefilepath)

end subroutine


    !------------------------------------------------------------------!
    !                     PLF CLOSED APPROXIMATION                     !
    !------------------------------------------------------------------!

    ! Initialization
    subroutine InitializePLF2DClosedApproximation(plf, ps)

        ! Description
        !============
        ! Initialize the general plf object. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedApproximationUDT)  :: plf 
        type(PolygonSetUDT), intent(in)             :: ps 

        ! Auxiliary
        class(PolygonLevelsetFunction2DUDT), allocatable :: plfe

        integer(I8)                         :: nx, ny
        real(R8), allocatable               :: xg(:), yg(:), vg(:), &
            vg2D(:, :), xgv(:), ygv(:)

        ! Initialize
        !===========
        ! Add polygonset structure
        plf%ps = ps 

        ! Construct closed exact levelset function
        call InitializePolygonLevelsetFunction2D(plfe, ps, &
            plf%options%optionsClosedExact)

        ! Construct 2D evaluation grid
        call plf%ConstructGrid(xg, yg, xgv, ygv, nx, ny)

        ! Evaluate exact representation
        allocate(vg(nx*ny))
        call plfe%Evaluate(xg, yg, 0, 0, vg)

        ! Reshape
        vg2D = reshape(vg, (/nx, ny/))

        ! Build polygon representation
        !=============================
        ! Construct
        call plf%interp%SetParameters(plf%options%meth, plf%options%C, plf%options%M)
        call plf%interp%Construct(xgv, ygv, vg2D)
        

    end subroutine

    ! Evaluation
    subroutine EvaluatePLF2DClosedApproximation(plf, xq, yq, derivx, &
        derivy, vq, varin, valuesin, dplfdvarin)

        ! Description
        !============
        ! Evaluate - simply call interpolant routine if no derivatives 
        ! need to be evaluate. Otherwise, if derivatives are needed, we
        ! apply the chain rule as follows:
        !
        !   dplf/dvar = dvq/db db/dvar
        !
        ! where dvq/db is the linearization of the interpolation routine
        ! w.r.t. the initial sample values 'b', and db/dvar the linearization
        ! of the initial sample values w.r.t. the desired variable. 
        ! Note that dvq/db is given by the interpolant as dvq/da da/db, 
        ! where da/db is the inverse of the system matrix used to 
        ! determine the interpolation weights a. db/dvar is the 
        ! linearization of the underlying closed exact polygon levelset
        ! function. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedApproximationUDT)  :: plf 
        real(R8), intent(in)            :: xq(:), yq(:)
        integer(I8), intent(in)         :: derivx, derivy
        real(R8), intent(out)           :: vq(size(xq))

        ! Optional arguments
        character(*), intent(in), optional      :: varin 
        real(R8), intent(in), optional          :: valuesin(:)
        type(MySparseUDT), optional             :: dplfdvarin

        character(:), allocatable               :: var 
        real(R8), allocatable                   :: values(:)
        type(MySparseUDT)                       :: dplfdvar

        ! Auxiliary
        class(PolygonLevelsetFunction2DUDT), allocatable :: plfe
        type(MySparseUDT)                       :: dinterpda, dadvqinit, &
             dvqinitdval, dvqdvqinit
        integer(I8)                             :: nx, ny
        real(R8), allocatable                   :: xg(:), &
            yg(:), vg(:), vg2D(:, :), xgv(:), ygv(:)
        logical                                 :: domemoryfriendly

        ! Check input
        !============
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
        if (present(dplfdvarin)) then 
            dplfdvar = dplfdvarin 
        end if 

        ! Set memory friendly option default to true
        domemoryfriendly = .true. 


        ! Call interpolant
        !=================
        ! Evaluate values
        call plf%interp%Evaluate(xq, yq, derivx, derivy, vq)


        select case (var)

        case ('no')

            ! No derivatives
            

        case ('polygonsetcoordinates')

            ! Only w.r.t. vessel coordinates. Need to differentiate 
            ! initialization and interpolation routines

            ! Construct closed exact levelset function
            call InitializePolygonLevelsetFunction2D(plfe, plf%ps, &
                plf%options%optionsClosedExact)

            ! Construct 2D evaluation grid
            call plf%ConstructGrid(xg, yg, xgv, ygv, nx, ny)

            ! Evaluate exact representation and its derivatives
            allocate(vg(size(xg, 1)*size(yg, 1)))
            call plfe%Evaluate(xg, yg, derivx, derivy, vg, 'polygonsetcoordinates', &
                values, dvqinitdval)

            ! Reshape
            vg2D = reshape(vg, (/nx, ny/))

            ! Compute derivative
            if (domemoryfriendly) then 
                ! Differentiate interpolated values w.r.t. input values (memory friendly)
                call plf%interp%EvaluateDiffInterp2Val(xq, yq, derivx, derivy, dvqdvqinit)

                ! Construct total derivative
                dplfdvar = dvqdvqinit*dvqinitdval

            else
                ! Memory unfriendly alternative
                ! Differentiate interpolant coefficient w.r.t. input values
                call plf%interp%EvaluateDiffCoef2Val(xgv, ygv, vg2D, dadvqinit)
                
                ! Construct
                call plf%interp%Construct(xgv, ygv, vg2D)

                ! Differentiate interpolant results w.r.t. coefficients
                call plf%interp%EvaluateDiffInterp2Coef(xq, yq, derivx, derivy, vq, &
                dinterpda)

                ! Construct total derivative
                dplfdvar = dinterpda*dadvqinit
                dplfdvar = dplfdvar*dvqinitdval
            end if 

        case default
            
            ! Throw error
            call gdErrorHandler('EvaluatePLF2DClosedApproximation: ' // & 
                ' variable not implemented')

        end select

        if (present(dplfdvarin)) then 
            dplfdvarin = dplfdvar 
        end if
 
    end subroutine

    ! Grid construction
    subroutine ConstructGridPLF2DClosedApproximation(plf, xg, yg, xgv, &
        ygv, nx, ny)

        ! Description
        !============
        ! Routine to create a 2D structured mesh for the closed 
        ! approximation plf function. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonLevelsetFunction2DClosedApproximationUDT)  :: plf 
        real(R8), intent(out), allocatable      :: xg(:), yg(:), &
            xgv(:), ygv(:)
        integer(I8), intent(out)                :: nx, ny

        ! Auxiliary
        real(R8)                            :: minx, maxx, miny, maxy, &
            dx, dy
        real(R8), allocatable               :: xp(:), yp(:)

        ! Loop
        integer(I8)                         :: k

        ! Initialize
        !===========
        ! Get polygon vertices 
        call plf%ps%GetVertices(xp, yp)

        ! Check interpolant range
        if (allocated(plf%options%xrange)) then 
            if (size(plf%options%xrange, 1) < 2) then 
                ! Set to default
                deallocate(plf%options%xrange)
                plf%options%xrange = xp
            end if 
            if (size(plf%options%yrange, 1) < 2) then 
                ! Set to default
                deallocate(plf%options%yrange)
                plf%options%yrange = yp
            end if 
        else
            ! Set to xp, yp
            plf%options%xrange = xp 
            plf%options%yrange = yp
        end if

        ! Evaluate exact representation
        !==============================
        ! Set dimensions of vertices (assumed resx, resy are number of 
        ! cells)
        nx = plf%options%resx + 1
        ny = plf%options%resy + 1

        ! Allocate
        allocate(xgv(nx), ygv(ny), xg(nx*ny), yg(nx*ny))

        ! Get polygon extent
        minx = minval(plf%options%xrange) 
        maxx = maxval(plf%options%xrange)
        miny = minval(plf%options%yrange) 
        maxy = maxval(plf%options%yrange)
        dx = maxx - minx 
        dy = maxy - miny

        ! Check how far the grid should extend
        minx = minx - plf%options%offsetx*dx 
        maxx = maxx + plf%options%offsetx*dx 
        miny = miny - plf%options%offsety*dy 
        maxy = maxy + plf%options%offsety*dy 
        
        ! Reset dx, dy
        dx = maxx - minx 
        dy = maxy - miny
        
        ! Gridding vectors
        xgv = [(k, k = 0, plf%options%resx)]*dx/plf%options%resx + minx 
        ygv = [(k, k = 0, plf%options%resy)]*dy/plf%options%resy + miny 
        
        call Construct2DStructuredGrid(xgv, ygv, nx, ny, xg, yg)

    end subroutine


end module