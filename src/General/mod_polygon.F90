!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module provides the 'polygonset' type (object), which can be used 
! to handle multiple polygons and polygon-like shapes. Each polygonset 
! contains a single or multiple polygons ('polygon' type)

! The coordinates may contain duplicate points, and the
! polygon(s) may intersect, be closed (this requires at least two duplicate
! points at the start and end), consist of multiple segments (separated by
! NaNs). Polygons can only be constructed by parsing x, y coordinate 
! lists, which may contain NaNs (IEEE standard) to distinguish 
! separate polygons. Due to the nature of the input,
! true branching polygons cannot be generated (which is often
! also not desired). However, by adding duplicate points, the polygon can
! in fact contain multiple times the same point and based on that self
! intersect and become non-simple. Note that a polygon is called 'simple'
! here if each vertex only occurs maximally twice in the polygon (so it can
! still intersect, simply not in a single node). 

! A polygonset object has the following fields
!       
!       np          number of distinct polygons
!       polygons    an np-by-1 structure containing np polygon objects
!       

! A polygon object has the following fields: 
!
!       nv, ne      number of unique polygon vertices (nv) and polygon 
!                   edges (ne)
!       nl          number of labels
!       x, y        vertex coordinates of the polygon (nv-by-1)
!       labels      labels for each (x, y) pair. Can be provided in the
!                   initialization (if not, default is nv-by-1 zero 
!                   array, but can be nv-by-nl)
!       edges       vertex pairs for each edge (ne-by-1, sorted in
!                   arbitrary direction)
!       vert        vertex ordering of the polygon (ne+1-by-1), may contain
!                   duplicates if the polygon is closed/non simple
!       ex, ey      polygon edge centers
!       tx, ty      polygon tangents (not normalized)
!       nx, ny      polygon normals (pointing to interior of polygon, not
!                   normalized)
!       tn, nn      tangent and normal length (for normalization if
!                   desired)
!       
!       isclosed    logical that indicates if the polygon is closed upon
!                   itself (i.e. the first and final edge have a vertex in
!                   common)
!       selfintersecting    logical that indicates if the polygon
!                           intersects with itself. Note: polygons that are
!                           closed but do not intersect anywhere else are
!                           not considered to be self intersecting. 
!       simplepolygon       logical that indicates if the polygon is a
!                           simple polygon, i.e. that each vertex belongs
!                           to maximally 2 edges. 

! The following routines are provided:
!       
!       ConstructPolygonSet(p, x, y) Polygon constructor, x and y should 
!                                   be 1D arrays of equal size, p is the
!                                   polygonset object 

! Notes
!======
! Note 1: a polygon should consist of at least one edge (i.e. 2 
! vertices). Otherwise, an error is thrown during construction. 

! Note 2: NaNs should be used to separate polygons. Therefore, they 
! should not occur in the first and last element of the coordinate 
! vectors (indicating an empty polygon, which is not supported, see 
! note 1). This is checked for and an error will be thrown. 

! Note 3: the intersection routines that are currently available can
! cope with most polygon behavior, except for exactly collinear edges.
! These will not be counted as intersections... 

module mod_polygon 

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_errorhandler
    use, intrinsic :: ieee_arithmetic
    use mod_plotter
    use mod_sort
    use mod_constants, only : pi_R8
    use omp_lib

    ! The usual
    implicit none
    save
    public 

    ! Constants
    real(R8), private, parameter        :: disttol = 1e-12 ! tolerance when computing distances
    real(R8), private, parameter        :: macheps = 1e-12 ! 'machine' precision

    !==================================================================!
    !                                                                  !
    !                                TYPES                             !
    !                                                                  !
    !==================================================================!

    ! The polygon type
    type PolygonUDT 

        ! The general polygon type, see documentation of this module for
        ! the fields and procedures. 

        integer(I8)                 :: nv, ne, nl

        logical                     :: isclosed, &
            issimple

        integer(I8), allocatable    :: edges(:, :), vert(:), labels(:, :)

        real(R8), allocatable       :: x(:), y(:), ex(:), ey(:), &
            tx(:), ty(:), tn(:), nx(:), ny(:), nn(:)
        

    contains 

        ! Construction 
        procedure :: Allocate       => AllocatePolygon
        procedure :: Deallocate     => DeallocatePolygon
        procedure, private  :: ConstructPolygon 
        procedure, private  :: ConstructPolygonNolabels
        generic   :: Construct      => ConstructPolygon, ConstructPolygonNolabels
        procedure :: Initialize     => InitializePolygon 

        ! Updating
        procedure :: UpdateCoordinates      => UpdatePolygonVertexCoordinates
        procedure :: Refine         => RefinePolygonUniform

        ! Metric computation & auxiliaries
        procedure :: ComputeMetrics => ComputePolygonMetrics
        procedure :: SetVert        => SetPolygonVertices
        procedure :: RemoveDuplicatePoints
        procedure :: IsClosedPolygon
        procedure :: IsSimplePolygon
        procedure :: IsSelfIntersectingPolygon
        procedure :: Inpolygon
        procedure :: Flip               => FlipPolygon
        procedure :: SelfIntersections  => PolygonSelfIntersections
        procedure :: GetSurfaceArea     => ComputePolygonSurfaceArea
        procedure, private  :: GetPolygonVertexID
        generic   :: GetVert            => GetPolygonVertexID

        ! Adaptations
        procedure :: Compress           => CompressPolygon  

    end type 

    ! The polygonset type
    type PolygonSetUDT 

        ! The general polygonset type, see documentation of this module
        ! for the fields and procedures

        integer(I8)                     :: np 
        type(PolygonUDT), allocatable   :: polygons(:)

    contains 

        ! Construction
        procedure, private  :: ConstructPolygonSetCoordinates
        procedure, private  :: ConstructPolygonSetFromPolygons
        procedure, private  :: ConstructPolygonSetFromEdges
        procedure, private  :: ConstructPolygonSetCoordinatesNolabels
        generic :: Construct        => ConstructPolygonSetCoordinates, &
            ConstructPolygonSetFromPolygons, ConstructPolygonSetFromEdges, &
            ConstructPolygonSetCoordinatesNolabels

        ! Operations
        procedure :: SelfIntersections      => PolygonSetSelfIntersections
        procedure :: OrientNestedClosedPolygons
        procedure :: Refine         => RefinePolygonSetUniform

        ! Data access
        procedure, private  :: GetPolygonSetEdgesCoordinates, &
            GetPolygonSetEdgesIDs
        generic   :: GetEdges               => GetPolygonSetEdgesCoordinates, &
            GetPolygonSetEdgesIDs
        procedure :: GetNormals             => GetPolygonSetNormals
        procedure :: GetTangents            => GetPolygonSetTangents
        procedure, private   :: GetPolygonSetVerticesCoordinates, &
            GetPolygonSetVerticesID
        generic   :: GetVertices            => GetPolygonSetVerticesCoordinates, &
            GetPolygonSetVerticesID
        procedure :: GetLabels              => GetPolygonSetVertexLabels
        procedure :: UpdateCoordinates      => UpdatePolygonSetVertexCoordinates

        ! Adaptations
        procedure :: Compress               => CompressPolygonSet          

        ! I/O
        procedure :: WriteData              => WritePolygonSetData        

    end type 

    !==================================================================!
    !                                                                  !
    !                             INTERFACES                           !
    !                                                                  !
    !==================================================================!

    ! Distance function for different inputs
    interface Distance 
        module procedure Distance0D
        module procedure Distance1D
    end interface


    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                         PolygonSet routines                      !
    !------------------------------------------------------------------!

    ! Construct the polygon set
    subroutine ConstructPolygonSetCoordinates(polygonset, x, y, labels)

        ! Description
        !============
        ! Construct a set of polygons from the given x and y 
        ! coordinates. NaNs should be present to indicate different 
        ! polygon pieces. First, the number of polygons are checked 
        ! by computing the number of NaNs etc. Afterwards, each polygon
        ! is constructed separately, including all metric and logical
        ! data. 

        ! Note: labels can be added (integer 2D array) to save other 
        ! data. These labels are propagated to the subsequent polygons
        ! and stored there locally. The labels do not need to contain
        ! NaNs at polygon boundaries. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(in)   :: x(:), y(:) 
        integer(I8), intent(in)             :: labels(:, :)
        class(PolygonSetUDT), intent(inout) :: polygonset 

        ! Auxiliary 
        integer(I8)                     :: nnans, nx, ny, nl1, nl2
        integer(I8), allocatable        :: nanloc(:), startind(:), &
            endind(:), nvpp(:) 

        real(R8), allocatable           :: tempx(:), tempy(:)
        integer(I8), allocatable        :: templabels(:, :)

        ! Loop
        integer(I8)                     :: i, k


        ! Initialize
        !===========
        ! Compute and check sizes
        nx = size(x)
        ny = size(y)
        nl1 = size(labels, 1)
        nl2 = size(labels, 2)
        if (nx .ne. ny) then 
            ! Call the error handler
            call PolygonErrorHandler('x and y should have equal number of elements')
        end if 
        if (nx .ne. nl1) then 
            ! Call error handles
            call PolygonErrorHandler('labels should have equal number of elements as coordinates')
        end if 

        ! Count NaNs and check
        nnans = count(isnan(x))
        allocate(nanloc(nnans))
        nanloc = pack([(k, k = 1, nx)], isnan(x))
        if (any(isnan(x) .neqv. isnan(y))) then 
            call PolygonErrorHandler('NaN values do not correspond in x and y coordinates')
        end if
        
        ! Determine polygons
        !===================
        ! Compute number of polygons
        polygonset%np = nnans + 1 

        ! Checks
        if (isnan(x(1)) .or. isnan(x(nx))) then  
            ! NaNs in first or last entry, throw error
            call PolygonErrorHandler('NaN values detected in first ' &
                // 'and/or last element of coordinate array. ' &
                // 'Not supported, check input')
        end if 

        ! Initialize polygons
        allocate(startind(polygonset%np), endind(polygonset%np), &
            nvpp(polygonset%np))
        startind    = [1, nanloc+1]
        endind      = [nanloc-1, nx] 
        nvpp        = endind - startind + 1
        if (allocated(polygonset%polygons)) then 
            deallocate(polygonset%polygons)
        end if 
        allocate(polygonset%polygons(polygonset%np))

        ! Check if all polygons have more than one vertex 
        if (any(nvpp < 2)) then 
            call PolygonErrorHandler('Polygons with less than two ' &
                // 'vertices found, not supported. Check input.')
        end if 
        
        ! Loop over all polygons
        do i = 1, polygonset%np 

            ! Build polygon
            tempx = x(startind(i):endind(i))
            tempy = y(startind(i):endind(i))
            templabels = labels(startind(i):endind(i), :)
            call polygonset%polygons(i)%Construct(tempx, tempy, templabels)

        end do 

        ! Housekeeping
        !=============
        deallocate(startind, endind, nvpp, nanloc)

    end subroutine

    ! Construct the polygon set (no labels)
    subroutine ConstructPolygonSetCoordinatesNolabels(polygonset, x, y)

        ! Description
        !============
        ! Construct a set of polygons from the given x and y 
        ! coordinates. NaNs should be present to indicate different 
        ! polygon pieces. First, the number of polygons are checked 
        ! by computing the number of NaNs etc. Afterwards, each polygon
        ! is constructed separately, including all metric and logical
        ! data. 

        ! Note: we simply call the routine with labels, but initialize
        ! the labels to an empty array

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(in)   :: x(:), y(:) 
        class(PolygonSetUDT), intent(inout)  :: polygonset 

        ! Auxiliary 
        integer(I8)                     :: nx
        integer(I8), allocatable        :: labels(:, :)

        ! Initialize
        !===========
        ! Set nx (dimension checks done in general routine)
        nx = size(x)

        ! Set labels
        allocate(labels(nx, 0))

        ! Call constructor
        !=================
        call polygonset%ConstructPolygonSetCoordinates(x, y, labels)

    end subroutine

    ! Construct starting from array of polygons
    subroutine ConstructPolygonSetFromPolygons(polygonset, polygons)

        ! Description
        !============
        ! Construct the set starting from already initialized polygons.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        type(PolygonUDT), intent(in)        :: polygons(:)
    
        ! Simply add polygons
        !====================
        ! Set number
        polygonset%np = size(polygons, 1)

        ! Assign
        polygonset%polygons = polygons


    end subroutine

    ! Construct starting from unsorted array of edges
    subroutine ConstructPolygonSetFromEdges(polygonset, edges, x, y, &
        extralabels_opt)

        ! Description
        !============
        ! Construct a set of polygons starting from an unsorted set of
        ! edges. These edges cannot contain NaNs. x, y are coordinate 
        ! vectors that should be indexable using the edge indices. 
        ! The vertex indices of the edges will be stored as labels in 
        ! the polygon structure.

        ! Note: optionally, additional vertex labels may be parsed, 
        ! which should also be the same (row) size as x

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)        :: polygonset 
        integer(I8), intent(in)     :: edges(:, :)
        real(R8), intent(in)        :: x(:), y(:)
        integer(I8), intent(in), optional   :: extralabels_opt(:, :)

        ! Auxiliary
        integer(I8)                 :: ne, np, tempne
        integer(I8), allocatable    :: sortindex(:), sortededges(:, :), &
            tempedges(:, :), tempvert(:), ps(:), pe(:), &
            templabels(:, :), extralabels(:, :)

        logical, allocatable        :: ispolygonstart(:), isbranchingpolygon(:)

        ! Loop
        integer(I8)                 :: i, j, k

        ! Initialize
        !===========
        ! Check edge dimensions
        if (size(edges, 2) /= 2) then 
            call gdErrorHandler('ConstructPolygonSetFromEdges: ' // &
                'edges should be ne-by-2 array')
        end if 

        ! Unpack
        ne = size(edges, 1)
        if (present(extralabels_opt)) then 
            extralabels = extralabels_opt 
        else
            allocate(extralabels(size(x), 0))
        end if 

        ! Allocate
        allocate(sortindex(ne), ispolygonstart(ne), isbranchingpolygon(ne))

        ! Extract polygon edges
        !======================
        ! Sort
        call SortPolygonEdges(edges, ne, sortindex, ispolygonstart, &
            isbranchingpolygon)
        allocate(sortededges, source=edges)
        sortededges = edges(sortindex, :)

        ! Get number of polygons
        np = count(ispolygonstart)

        ! Construct polygon set
        !======================
        ! Allocate
        polygonset%np = np 
        allocate(polygonset%polygons(np))

        ! Construct polygon start and end indices
        allocate(ps(np), pe(np))
        ps = pack([(k, k = 1, ne)], ispolygonstart)
        pe = [ps(2:np)-1, ne]

        ! Construct polygons
        do i = 1, np
            ! Get polygon edges
            tempne = pe(i) - ps(i) + 1
            tempedges = sortededges(ps(i):pe(i), :)

            ! Get polygon vertices
            allocate(tempvert(tempne+1))
            call ExtractPolygonVertices(tempedges, tempne, &
                tempvert)

            ! Set labels as vertex indices
            allocate(templabels(tempne+1, 1 + size(extralabels, 2)))
            templabels(:, 1) = tempvert
            do j = 1, size(extralabels, 2)
                templabels(:, j+1) = extralabels(tempvert, j)
            end do 

            ! Construct polygon
            call polygonset%polygons(i)%Construct(x(tempvert), &
                y(tempvert), templabels)

            ! Housekeeping
            deallocate(tempvert, templabels)
        end do

        

    end subroutine

    ! Polygon set refiner (uniform)
    subroutine RefinePolygonSetUniform(polygonset, dlmax, minsplit)

        ! Description
        !============
        ! This routine refines each of the polygons in the polygonset 
        ! in a uniform way. This routine is just a wrapper for the 
        ! polygon refinement routine - see that routine for additional
        ! information. 


        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        real(R8), intent(in)                :: dlmax 
        integer(I8), intent(in)             :: minsplit

        ! Auxiliary

        ! Loop
        integer(I8)                         :: i 

        ! Refine
        !=======
        ! Simply call polygon refiner for each polygon...
        do i = 1, polygonset%np
            call polygonset%polygons(i)%Refine(dlmax, minsplit)
        end do

    end subroutine

    ! Polygon set compressor 
    subroutine CompressPolygonSet(polygonset)

        ! Description
        !============
        ! 'Compress' the polygon set by removing points that do not 
        ! contribute to the shape of the polygon. Primarily useful to 
        ! simplify polygons that were refined before. Differences in the
        ! order of machine precision may originate from using this 
        ! function. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset

        ! Loop
        integer(I8)                         :: i 

        ! Call polygon compressor
        !========================
        do i = 1, polygonset%np
            call polygonset%polygons(i)%Compress()
        end do 

    end subroutine

    ! Update polygonset coordinates
    subroutine UpdatePolygonSetVertexCoordinates(polygonset, xp, yp)

        ! Description
        !============
        ! Update the vertex coordinates of a polygonset. First, it is 
        ! checked that xp, yp have the correct dimensions (i.e. equal 
        ! to the number of vertices of all polygons). Then, we loop over
        ! the polygons and attribute the correct vertices. It is assumed
        ! that the polygons are correctly updated by attributing these
        ! vertices with the polygon specific subroutine. Note that 
        ! self-intersections or polygon closedness may change!

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)            :: polygonset 
        real(R8), intent(in)            :: xp(:), yp(:)

        ! Auxiliary
        integer(I8)                     :: nvtot
        integer(I8), allocatable        :: vIDs(:), tempind(:)

        ! Loop
        integer(I8)                     :: i, j, nv

        ! Initialize
        !===========
        ! Associate
        associate(&
            np          => polygonset%np,       &
            pol         => polygonset%polygons  &
            )

        ! Compute total number of vertices
        call polygonset%GetVertices(vIDs)
        nvtot = size(vIDs, 1)

        ! Check
        if ((size(xp, 1) /= nvtot) .or. (size(xp, 1) /= size(yp, 1))) then 
            ! Throw error
            call gdErrorHandler('UpdatePolygonSetVertexCoordinates: ' // &
                'xp, yp do not have as many elements as there are vertices in the ' // &
                'polygon set. Check input.')
        end if 

        ! Update coordinates
        !===================
        nv = 0
        do i = 1, np
            ! Set index vector
            tempind = [(j, j = nv+1, nv+pol(i)%nv)]

            ! Update
            call pol(i)%UpdateCoordinates(xp(tempind), yp(tempind))

            ! Update counter
            nv = nv + pol(i)%nv 
        end do
        
        ! Housekeeping
        !=============
        end associate
        

    end subroutine

    ! Compute polygon set self-intersections
    subroutine PolygonSetSelfIntersections(polygonset, x, y, p1, p2, &
        s1, s2)

        ! Description
        !============
        ! Compute the self-intersections of a polygon set. Here, both
        ! theintersections between different polygons of the set
        ! and self-intersections of polygons are computed. 
        ! As output, the routine
        ! returns the x, y coordinates of the intersections, the indices
        ! of the polygons of these intersections, and the indices of the
        ! segments in said polygons of these intersections. The routine
        ! builds upon the PolygonIntersections routine for computing
        ! the intersections of two polygons, and the 
        ! PolygonSelfIntersections routine for computing 
        ! self-intersections of a single polygon

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                    :: polygonset 
        real(R8), allocatable, intent(out)      :: x(:), y(:)
        integer(I8), allocatable, intent(out)   :: p1(:), p2(:), &
            s1(:), s2(:)
        
        ! Auxiliary
        integer(I8)                             :: ni, counter, sz, &
            szold, szmult 
        real(R8), allocatable                   :: tempx(:), tempy(:), &
            xi(:), yi(:) 
        integer(I8), allocatable                :: temps1(:), &
            temps2(:), si1(:), si2(:), tempp1(:), tempp2(:)

        ! Loop
        integer(I8)                             :: i, j 

        ! Memory mgmt
        integer(I8), allocatable                :: mgmti(:, :)
        real(R8), allocatable                   :: mgmtr(:, :)

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if
        if (allocated(s1)) then 
            deallocate(s1) 
        end if
        if (allocated(s2)) then 
            deallocate(s2) 
        end if
        if (allocated(p1)) then 
            deallocate(p1) 
        end if
        if (allocated(p2)) then 
            deallocate(p2) 
        end if

        ! Associate
        associate(&
            np          => polygonset%np,       &
            p           => polygonset%polygons)

        ! Initialize
        counter     = 0 ! intersection counter 
        ni          = 0
        szold       = 0
        sz          = 2 ! initial size of intersection array
        szmult      = 2 ! size multiplier 

        ! Allocate
        allocate(tempx(sz), tempy(sz), temps1(sz), temps2(sz), &
            tempp1(sz), tempp2(sz))

        ! Compute intersections between polygons
        !=======================================
        do i = 1, np-1  ! Loop over all polygons-1
            do j = i+1, np ! Loop over remaining polygons
                ! Compute intersections
                call PolygonIntersections(p(i), p(j), xi, yi, si1, si2)

                ! Check if intersections were found
                ni = size(xi)
                if (ni > 0) then 
                    ! Memory MGMT
                    if (counter + ni > sz) then 
                        ! Store old size
                        szold = sz

                        ! Store old values
                        allocate(mgmti(szold, 4), mgmtr(szold, 2))
                        mgmti(:, 1) = temps1 
                        mgmti(:, 2) = temps2   
                        mgmti(:, 3) = tempp1 
                        mgmti(:, 4) = tempp2
                        mgmtr(:, 1) = tempx 
                        mgmtr(:, 2) = tempy

                        ! Adjust size
                        do while (sz < counter+ni)
                            sz = sz*szmult 
                        end do

                        ! Reallocate
                        deallocate(tempx, tempy, temps1, temps2, &
                            tempp1, tempp2) 
                        allocate(tempx(sz), tempy(sz), temps1(sz), &
                            temps2(sz), tempp1(sz), tempp2(sz))

                        ! Add
                        tempx(1:szold) = mgmtr(:, 1) 
                        tempy(1:szold) = mgmtr(:, 2)
                        temps1(1:szold) = mgmti(:, 1)
                        temps2(1:szold) = mgmti(:, 2)
                        tempp1(1:szold) = mgmti(:, 3)
                        tempp2(1:szold) = mgmti(:, 4)

                        ! Deallocate mgmt arrays
                        deallocate(mgmti, mgmtr)
                    end if 

                    ! Add intersections
                    tempx(counter+1:counter+ni) = xi 
                    tempy(counter+1:counter+ni) = yi
                    temps1(counter+1:counter+ni) = si1 
                    temps2(counter+1:counter+ni) = si2
                    tempp1(counter+1:counter+ni) = i 
                    tempp2(counter+1:counter+ni) = j

                    ! Update counter
                    counter = counter + ni
                end if 
            end do  
        end do

        ! Compute polgon self-intersections
        !==================================
        do i = 1, np 
            ! Compute intersections
            call PolygonSelfIntersections(p(i), xi, yi, si1, si2)

            ! Check if intersections were found
            ni = size(xi)
            if (ni > 0) then 
                ! Memory MGMT
                if (counter + ni > sz) then 
                    ! Store old size
                    szold = sz

                    ! Store old values
                    allocate(mgmti(szold, 4), mgmtr(szold, 2))
                    mgmti(:, 1) = temps1 
                    mgmti(:, 2) = temps2   
                    mgmti(:, 3) = tempp1 
                    mgmti(:, 4) = tempp2
                    mgmtr(:, 1) = tempx 
                    mgmtr(:, 2) = tempy

                    ! Adjust size
                    do while (sz < counter+ni)
                        sz = sz*szmult 
                    end do

                    ! Reallocate
                    deallocate(tempx, tempy, temps1, temps2, &
                        tempp1, tempp2) 
                    allocate(tempx(sz), tempy(sz), temps1(sz), &
                        temps2(sz), tempp1(sz), tempp2(sz))

                    ! Add
                    tempx(1:szold) = mgmtr(:, 1) 
                    tempy(1:szold) = mgmtr(:, 2)
                    temps1(1:szold) = mgmti(:, 1)
                    temps2(1:szold) = mgmti(:, 2)
                    tempp1(1:szold) = mgmti(:, 3)
                    tempp2(1:szold) = mgmti(:, 4)

                    ! Deallocate mgmt arrays
                    deallocate(mgmti, mgmtr)
                end if 

                ! Add intersections
                tempx(counter+1:counter+ni) = xi 
                tempy(counter+1:counter+ni) = yi
                temps1(counter+1:counter+ni) = si1 
                temps2(counter+1:counter+ni) = si2
                tempp1(counter+1:counter+ni) = i 
                tempp2(counter+1:counter+ni) = i

                ! Update counter
                counter = counter + ni
            end if 
        end do  

        ! Add to output
        allocate(x(counter), y(counter), s1(counter), s2(counter), &
            p1(counter), p2(counter)) 
        x = tempx(1:counter) 
        y = tempy(1:counter) 
        s1 = temps1(1:counter) 
        s2 = temps2(1:counter) 
        p1 = tempp1(1:counter)
        p2 = tempp2(1:counter)

        ! Housekeeping
        !=============
        end associate 

        deallocate(tempx, tempy, temps1, temps2, tempp1, tempp2)



    end subroutine

    ! Orient nested closed polygons
    subroutine OrientNestedClosedPolygons(polygonset, flag, polygonlevels)

        ! Description
        !============
        ! This routine checks first whether the polygon set consists of 
        ! closed, (non-intersecting) nested polygons. If this is the 
        ! case, the routine orients them as such that the interior of 
        ! the domain is well defined.

        ! This routine orients nested polygons such that the interior of the domain
        ! is well defined. The initial orientation is such that the surface area
        ! computed by integrating each line segment over the x-axis is positive.
        ! Each nested polygon is then reoriented by changing the vertex order of
        ! that polygon. 

        ! Algorithm
        !==========
        ! To sort the polygons, we first determine which polygons are contained by
        ! others and store this information in a connectivity matrix
        ! 'inpolygonmatrix'. If the (i, j)th element is nonzero, this means that
        ! the j-th polygon lies within the i-th polygon. Afterwards, we have the
        ! following algorithm to determine the orientation:
        !
        ! 1) determine the initial orientation (this is taken equal to the sign of
        ! the signed surface area of the polygon).
        ! 2) For each polygon that does not have any parents, check if the
        ! orientation matches the desired orientation. If not, change the sequence
        ! of coordinates and related quantities (TPind etc)
        ! 3) Mark these polygons as sorted and remove them as parents from any
        ! other polygon (this can be done by modifying the inpolygonmatrix). Change
        ! the sign of the desired orientation (this ensures a proper definition of
        ! the interior of the domain)
        ! 4) Repeat 2-3 until all polygons are found

        ! Notes
        !======
        ! Note 1: this routine has not yet been thoroughly verified for 
        ! multiple polygons - a warning is thrown 

        ! Note 2: if the routine exits successfully, the output flag 
        ! is equal to zero. If an exception occurs, it is larger than
        ! zero, and typically a warning message will be printed. Some
        ! flag meanings:
        ! 0:    success
        ! 1:    exit due to no polygons present in the set
        ! 2:    intersecting closed polygons
        ! 3:    coinciding closed polygons  
        ! 4:    open polygons present
        ! 5:    self-intersecting polygons present

        ! Note 3: if desired, also the 'level' of the polygon is returned,
        ! where the level indicates which nestedness level the polygon 
        ! has. Here, level 0 indicates polygons that do not lie in any
        ! other polygon, level 1 indicates polygons that lie in one 
        ! upper polygon, etc. If for any reason the routine does not 
        ! succeed, all levels are set to -1. Note that the polygon level
        ! is easily obtained from the inpolygonmatrix, as this is simply
        ! the column sum of that matrix (where true -> 1, false -> 0)

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)        :: polygonset 
        integer(I8)                 :: flag
        integer(I8), allocatable, dimension(:), intent(out), optional :: &
            polygonlevels

        ! Auxiliary
        integer(I8)                 :: orientation, nv

        real(R8), allocatable       :: polygonarea(:), ipx(:), ipy(:), &
            yf(:), dx(:), jpx(:), jpy(:)   

        logical, allocatable        :: inpolygonmatrix(:, :), jini(:), &
            iinj(:), ispolygonfound(:), doflip(:), dopolyg(:)

        ! Loop
        integer(I8)                 :: i, j 

        ! Initialize
        !===========
        ! Set flag (zero for success, > 0 for failure)
        flag = 1

        ! Set optional output if present
        if (present(polygonlevels)) then 
            if (allocated(polygonlevels)) deallocate(polygonlevels)
            allocate(polygonlevels(polygonset%np))
            polygonlevels = -1 ! initial value
        end if 

        ! Associate
        associate( &
                np  => polygonset%np, &
                p   => polygonset%polygons)

        ! Check polygon status
        do i = 1, np
            ! Check closedness
            if (.not. p(i)%isclosed) then 
                ! Adjust flag and exit
                call PolygonWarningHandler(&
                    'OrientNestedClosedPolygons: open polygons detected. Returning...')
                flag = 4
                return 
            end if 
            if (p(i)%IsSelfIntersectingPolygon()) then 
                ! Adjust flag and exit
                call PolygonWarningHandler(&
                    'OrientNestedClosedPolygons: self-intersecting polygons detected. Returning...')
                flag = 5
                return 
            end if 
        end do  

        ! Exit in trivial cases
        if (np == 0) then 
            return 
        end if 

        ! Determine polygon nestedness
        !=============================
        ! Initialize
        allocate(polygonarea(np), inpolygonmatrix(np, np))
        polygonarea(:) = 0
        inpolygonmatrix = .false.

        ! Compute area and nestedness
        do i = 1, np
            ! Initialize
            associate(&
                ne      => p(i)%ne)
            nv = ne+1

            allocate(ipx(nv), ipy(nv), yf(ne), dx(ne))

            ! Get current polygon coordinates
            ipx = p(i)%x(p(i)%vert) 
            ipy = p(i)%y(p(i)%vert) 

            ! Compute the surface area (basically integral over x with midpoint
            ! rule)
            yf = 0.5*(p(i)%y(p(i)%vert(1:ne)) + p(i)%y(p(i)%vert(2:ne+1)))
            dx = (p(i)%x(p(i)%vert(2:ne+1)) - p(i)%x(p(i)%vert(1:ne)))
            polygonarea(i) = sum(yf*dx)

            ! Determine which polygons lie within this polygon or in which polygons
            ! this polygon lies
            do j = 1, np
                ! Skip for same polygons
                if (j == i) then 
                    ! Skip
                    cycle
                end if

                ! Initialize
                associate(&
                    nvj      => p(j)%nv)
                allocate(jpx(nvj), jpy(nvj))

                ! Get next polygon coordinates
                jpx = p(j)%x
                jpy = p(j)%y

                ! Check if in polygon
                call p(i)%inpolygon(jpx, jpy, jini)
                call p(j)%inpolygon(ipx, ipy, iinj)

                ! Sanity checks
                if ((.not. all(jini) .and. (.not. all(.not. jini))) .or. &
                     (.not. all(iinj) .and. (.not. all(.not. iinj))) ) then 
                    ! Some polygons are only partly in another polygon. This would
                    ! indicate intersections between closed polygons, which is not
                    ! supported (and should in fact already be flagged above).
                    ! Call warning handler and exit
                    call PolygonWarningHandler(&
                        'OrientNestedClosedPolygons: when sorting ' // &
                        ' polygons, some polygons seem to ' // &
                        'intersect. Returning...')
                    flag = 2
                    return 
                end if
                if (all(jini) .and. all(iinj)) then 
                    ! Coinciding polygons - also not supported
                    call PolygonWarningHandler(&
                        'OrientNestedClosedPolygons: when sorting ' // &
                        'polygons, detected coinciding polygons, ' // &
                        'not supported. Returning...')
                    flag = 3
                    return 
                end if

                ! Add to matrix
                if (all(jini)) then 
                    inpolygonmatrix(i, j) = .true.
                elseif (all(iinj)) then 
                    inpolygonmatrix(j, i) = .true.
                end if 

                !    Housekeeping
                end associate 
                deallocate(jpx, jpy)
            end do 

            !    Housekeeping
            end associate 
            deallocate(ipx, ipy, yf, dx)
        end do

        ! Set levels
        !===========
        if (present(polygonlevels)) then 
            do i = 1, np
                polygonlevels(i) = count(inpolygonmatrix(:, i))
            end do 
        end if 
        
        ! Check orientation
        !==================
        ! Initialize
        allocate(ispolygonfound(np), doflip(np), dopolyg(np))
        ispolygonfound(:)   = .false.
        orientation         = 1
        doflip(:)           = .false.

        ! Loop
        do while (.not. all(ispolygonfound)) 
            ! Find the next overarching polygon
            dopolyg(:) = .false. 
            do i = 1, np
                if (.not. ispolygonfound(i)) then 
                    if (.not. any(inpolygonmatrix(:, i))) then  
                        ! Polygon i does not lie in any other polygon
                        dopolyg(i) = .true.
                        ispolygonfound(i) = .true.
                    end if
                end if
            end do
            
            ! Remove
            do i = 1, np 
                if (ispolygonfound(i)) then 
                    inpolygonmatrix(i, :) = .false. 
                end if 
            end do 

            ! Check orientation
            do i = 1, np
                if (dopolyg(i)) then 
                    ! Check orientation
                    if (polygonarea(i)*orientation < 0) then 
                        ! If the signs differ, the expression is 
                        ! negative and therefore the polygon should
                        ! be flipped
                        doflip(i) = .true.
                    end if
                end if
            end do

            ! Adjust orientation
            orientation = -1*orientation;
        end do 

        ! Flip the polygon fields - sadly, no dynamic field naming is 
        ! possible for UDTs in fortran AFAIK...
        do i = 1, np
            if (doflip(i)) then 
                ! Flip
                call p(i)%Flip()
            end if 
        end do 

        ! If we got here, everything worked fine. Set flag to zero
        flag = 0

        ! Housekeeping
        !=============
        end associate

        



    end subroutine

    ! Get polygonset edges
    subroutine GetPolygonSetEdgesCoordinates(polygonset, xe, ye)

        ! Description
        !============
        ! Return the edge coordinates of the entire polygon set - useful
        ! for bulk geometric operations such as normal computations etc.
        ! It is assumed that the polygonset is fully up to date. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        real(R8), allocatable, intent(out)  :: xe(:, :), ye(:, :)

        ! Auxiliary
        integer(I8)                         :: ne, ce  

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Deallocate if already allocated
        if (allocated(xe)) then 
            deallocate(xe)
        end if 
        if (allocated(ye)) then 
            deallocate(ye) 
        end if 

        ! Build edges
        !============
        ! Associate
        associate( &
            pol         => polygonset%polygons)
        
        ! Precompute the total number of edges
        ne = 0
        do i = 1, polygonset%np 
            ne = ne + polygonset%polygons(i)%ne             
        end do 

        ! Allocate
        allocate(xe(ne, 2), ye(ne, 2))

        ! Loop and add
        ce = 0 ! edge counter
        do i = 1, polygonset%np 
            ! Add coordinates
            xe(ce+1:ce+pol(i)%ne, 1) = pol(i)%x(pol(i)%edges(:, 1))
            xe(ce+1:ce+pol(i)%ne, 2) = pol(i)%x(pol(i)%edges(:, 2))
            ye(ce+1:ce+pol(i)%ne, 1) = pol(i)%y(pol(i)%edges(:, 1))
            ye(ce+1:ce+pol(i)%ne, 2) = pol(i)%y(pol(i)%edges(:, 2))

            ! Update counter
            ce = ce + pol(i)%ne 
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    subroutine GetPolygonSetEdgesIDs(polygonset, edges, polygonID)

        ! Description
        !============
        ! Similar to GetPolygonSetEdgesCoordinates, but now returns the
        ! edges by their vertex ID pairs. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                    :: polygonset 
        integer(I8), allocatable, intent(out)   :: edges(:, :)
        integer(I8), allocatable, intent(out), optional :: polygonID(:)

        ! Auxiliary
        integer(I8)                             :: netot 
        integer(I8), allocatable, dimension(:)  :: pID 

        ! Loop
        integer(I8)                             :: i, vc, ec

        ! Initialize
        !===========
        associate(&
            np          => polygonset%np,       &
            pol         => polygonset%polygons  &
            )

        ! Compute
        !========
        ! Precompute sizes
        netot = 0
        do i = 1, np 
            netot = netot + pol(i)%ne 
        end do

        ! Compute edges
        allocate(edges(netot, 2), pID(netot))
        vc = 0 ! keep track of vertices to update local polygon vertex IDs
        ec = 0
        do i = 1, np 
            ! Add
            edges(ec+1:ec+pol(i)%ne, :) = pol(i)%edges + vc 
            pID(ec+1:ec+pol(i)%ne)  = i 

            ! Update counters
            vc = vc + pol(i)%nv 
            ec = ec + pol(i)%ne
        end do 

        ! Housekeeping
        !=============
        if (present(polygonID)) then 
            polygonID = pID 
        end if 
        end associate
        

    end subroutine

    ! Get polygonset tangents (not normalized)
    subroutine GetPolygonSetTangents(polygonset, tx, ty, tn)

        ! Description
        !============
        ! Return the edge tangents of the entire polygon set - useful
        ! for bulk geometric operations such as normal computations etc.
        ! It is assumed that the polygonset is fully up to date. 
        ! Note that the tangents are not normalized, hence we return 
        ! the tangent length in tn 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        real(R8), allocatable, intent(out)  :: tx(:), ty(:), tn(:)

        ! Auxiliary
        integer(I8)                         :: ne, ce  

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Deallocate if already allocated
        if (allocated(tx)) then 
            deallocate(tx)
        end if 
        if (allocated(ty)) then 
            deallocate(ty) 
        end if 
        if (allocated(tn)) then 
            deallocate(tn) 
        end if

        ! Build edges
        !============
        ! Associate
        associate( &
            pol         => polygonset%polygons)
        
        ! Precompute the total number of edges
        ne = 0
        do i = 1, polygonset%np 
            ne = ne + polygonset%polygons(i)%ne             
        end do 

        ! Allocate
        allocate(tx(ne), ty(ne), tn(ne))

        ! Loop and add
        ce = 0 ! edge counter
        do i = 1, polygonset%np 
            ! Add coordinates
            tx(ce+1:ce+pol(i)%ne)   = pol(i)%tx 
            ty(ce+1:ce+pol(i)%ne)   = pol(i)%ty 
            tn(ce+1:ce+pol(i)%ne)   = pol(i)%tn 

            ! Update counter
            ce = ce + pol(i)%ne 
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Get polygonset normals (not normalized)
    subroutine GetPolygonSetNormals(polygonset, nx, ny, nn)

        ! Description
        !============
        ! Return the edge normals of the entire polygon set - useful
        ! for bulk geometric operations such as normal computations etc.
        ! It is assumed that the polygonset is fully up to date. 
        ! Note that the normals are not normalized, hence we return 
        ! the normal length in nn 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        real(R8), allocatable, intent(out)  :: nx(:), ny(:), nn(:)

        ! Auxiliary
        integer(I8)                         :: ne, ce  

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Deallocate if already allocated
        if (allocated(nx)) then 
            deallocate(nx)
        end if 
        if (allocated(ny)) then 
            deallocate(ny) 
        end if 
        if (allocated(nn)) then 
            deallocate(nn) 
        end if

        ! Build edges
        !============
        ! Associate
        associate( &
            pol         => polygonset%polygons)
        
        ! Precompute the total number of edges
        ne = 0
        do i = 1, polygonset%np 
            ne = ne + polygonset%polygons(i)%ne             
        end do 

        ! Allocate
        allocate(nx(ne), ny(ne), nn(ne))

        ! Loop and add
        ce = 0 ! edge counter
        do i = 1, polygonset%np 
            ! Add coordinates
            nx(ce+1:ce+pol(i)%ne)   = pol(i)%nx 
            ny(ce+1:ce+pol(i)%ne)   = pol(i)%ny 
            nn(ce+1:ce+pol(i)%ne)   = pol(i)%nn 

            ! Update counter
            ce = ce + pol(i)%ne 
        end do 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Get polygonset Vertices 
    subroutine GetPolygonSetVerticesCoordinates(polygonset, xp, yp, polygonID)

        ! Description
        !============
        ! Return all coordinates of the Vertices in an unspecified order.
        ! The points should in principle be unique. Useful for 
        ! operations that only considers the points and not the edges
        ! of the polygon in any particular order. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                :: polygonset 
        real(R8), allocatable, intent(out)  :: xp(:), yp(:)
        integer(I8), allocatable, intent(out), optional     :: polygonID(:)

        ! Auxiliary
        integer(I8)                         :: nv
        integer(I8), allocatable, dimension(:)  :: pID

        ! Loop
        integer(I8)                         :: i, ce

        ! Initialize
        !===========
        ! Deallocate if already allocated
        if (allocated(xp)) then 
            deallocate(xp)
        end if 
        if (allocated(yp)) then 
            deallocate(yp) 
        end if 

        ! Build points
        !=============
        ! Associate
        associate( &
            pol         => polygonset%polygons)
        
        ! Precompute the total number of points
        nv = 0
        do i = 1, polygonset%np 
            nv = nv + polygonset%polygons(i)%nv             
        end do 

        ! Allocate
        allocate(xp(nv), yp(nv), pID(nv))

        ! Loop and add
        ce = 0 ! edge counter
        do i = 1, polygonset%np 
            ! Add coordinates
            xp(ce+1:ce+pol(i)%nv)   = pol(i)%x 
            yp(ce+1:ce+pol(i)%nv)   = pol(i)%y 
            pID(ce+1:ce+pol(i)%nv)  = i

            ! Update counter
            ce = ce + pol(i)%nv 
        end do 

        ! Housekeeping
        !=============
        if (present(polygonID)) then 
            polygonID = pID 
        end if 
        end associate

    end subroutine

    subroutine GetPolygonSetVerticesID(polygonset, ID, polygonID)

        ! Description
        !============
        ! Return the ID vector of all points in the polygon set (IDs). 
        ! Additionally, the polygon index may be returned as optional
        ! output argument. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                    :: polygonset 
        integer(I8), allocatable, intent(out)   :: ID(:)
        integer(I8), allocatable, intent(out), optional     :: polygonID(:) 

        ! Auxiliary
        integer(I8)                         :: nv  
        integer(I8), allocatable            :: tID(:), pID(:)

        ! Loop
        integer(I8)                         :: i, vc

        ! Initialize
        !===========
        ! Deallocate if already allocated
        if (allocated(ID)) then 
            deallocate(ID)
        end if 

        ! Build points
        !=============
        ! Associate
        associate( &
            pol         => polygonset%polygons)
        
        ! Precompute the total number of points
        nv = 0
        do i = 1, polygonset%np 
            nv = nv + polygonset%polygons(i)%nv             
        end do 

        ! Allocate
        allocate(ID(nv), pID(nv))

        ! Loop and add
        vc = 0 ! vertex counter
        do i = 1, polygonset%np 
            ! Add coordinates
            call pol(i)%GetVert(tID)
            ID(vc+1:vc+pol(i)%nv)   = tID + vc ! update
            pID(vc+1:vc+pol(i)%nv)  =  i 

            ! Update counter
            vc = vc + pol(i)%nv 
        end do 

        ! Housekeeping
        !=============
        if (present(polygonID)) then 
            polygonID = pID 
        end if 
        end associate

    end subroutine

    ! Get polygonset vertex labels
    subroutine GetPolygonSetVertexLabels(polygonset, labels)

        ! Description
        !============
        ! Return the polygonset vertex labels in the same sequence
        ! as the vertex coordinates. Note that in general the number of
        ! labels between polygons doesn't have to be the same. If this 
        ! is the case, the labels that are returned are sized to the
        ! maximum number of labels and non-used labels are filled with
        ! zeros.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)            :: polygonset 
        integer(I8), allocatable        :: labels(:, :)

        ! Auxiliary
        integer(I8)                     :: maxnl, nv 

        ! Loop
        integer(I8)                     :: i, lc

        ! Initialize
        !===========
        ! Associate
        associate(pol   => polygonset%polygons)

        ! Compute maximal number of labels and number of vertices
        maxnl = 0
        nv = 0
        do i = 1, polygonset%np 
            maxnl = max(maxnl, size(pol(i)%labels, 2))
            nv = nv + pol(i)%nv
        end do 

        ! Allocate
        allocate(labels(nv, maxnl))

        ! Extract labels
        !===============
        lc = 0
        do i = 1, polygonset%np 
            ! Assign
            labels(lc+1:lc+pol(i)%nv, 1:size(pol(i)%labels, 2)) = &
                pol(i)%labels

            ! Update counter
            lc = lc + pol(i)%nv 
        end do

        ! Housekeeping
        !=============
        end associate 

    end subroutine

    ! Write polygonset vertex data
    subroutine WritePolygonSetData(polygonset, filepath)

        ! Description
        !============
        ! This routine writes out the polygon vertex data in x, y format
        ! for all polygons it contains. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonSetUDT)                    :: polygonset 
        character(*),  intent(in)               :: filepath

        ! Auxiliary
        integer(I8)                             :: ndata 
        real(R8)                                :: nan 
        real(R8), allocatable                   :: tempx(:), tempy(:)

        ! Loop
        integer(I8)                             :: i, cc

        ! Extract polygon data
        !=====================
        ! Unpack
        associate( &
            np          => polygonset%np,       &
            pol         => polygonset%polygons)

        ! Compute number of data entries
        ndata = 0
        do i = 1, np 
            ! Vertices of polygon
            ndata = ndata + pol(i)%ne+1
        end do 
        
        ! Account for NaNs
        ndata = ndata + np - 1

        ! Allocate
        allocate(tempx(ndata), tempy(ndata))

        ! Loop
        cc = 0
        do i = 1, np 
            ! Add coordinates
            tempx(cc+1:cc+pol(i)%ne+1) = pol(i)%x(pol(i)%vert)
            tempy(cc+1:cc+pol(i)%ne+1) = pol(i)%y(pol(i)%vert)

            ! Update counter
            cc = cc + pol(i)%ne + 1

            ! Add NaN
            if (i < np) then 
                tempx(cc+1) = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                tempy(cc+1) = IEEE_VALUE(nan, IEEE_QUIET_NAN)

                ! Update counter
                cc = cc + 1
            end if 

        end do

        ! Write
        !======
        call Write2DPolygonData(tempx, tempy, filepath)

        ! Housekeeping
        !=============
        end associate 
        deallocate(tempx, tempy)
        
        
    end subroutine

    !------------------------------------------------------------------!
    !                    PolygonSet derivative routines                !
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    !                          Polygon routines                        !
    !------------------------------------------------------------------!

    ! Constructor
    subroutine ConstructPolygon(polygon, x, y, labels, forcevert_opt) 

        ! Description
        !============
        ! Construct a single polygon from the given x and y coordinates.
        ! It is assumed that there are no NaNs present anymore in the 
        ! coordinates and that the sizes of x and y are the same. 

        ! Note: to prevent the duplicate point remover to remove certain
        ! key points, one can specify those not to be removed in 
        ! forcevert_opt (this should be a logical). By default, any point
        ! can be removed. This logical vector should have the same number
        ! of elements as the initial x, y coordinates

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                       :: polygon 
        real(R8), intent(in)                    :: x(:), y(:)
        integer(I8), intent(in)                 :: labels(:, :) 
        logical, intent(in), optional           :: forcevert_opt(:)

        ! Auxiliary

        ! Loop

        ! Construct
        !==========
        ! Set coordinates & initialize vertices etc
        call polygon%Initialize(x, y, labels)

        ! Check for duplicate points, remove
        call polygon%RemoveDuplicatePoints(forcevert_opt)

        ! Extract the vertices
        call polygon%SetVert()

        ! Check if the polygon is closed
        call polygon%IsClosedPolygon()

        ! Construct the polygon metrics
        call polygon%ComputeMetrics()

        ! Check if the polygon is simple
        call polygon%IsSimplePolygon()

        ! Check if the polygon self intersects
        ! call polygon%IsSelfIntersectingPolygon()


    end subroutine

    ! Constructor, without labels
    subroutine ConstructPolygonNolabels(polygon, x, y, forcevert_opt)

        ! Description
        !============
        ! Wrapper for construction without labels. Labels are simply
        ! initialized to a zero array 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)           :: polygon 
        real(R8), intent(in)        :: x(:), y(:)
        logical, intent(in), optional           :: forcevert_opt(:)

        ! Auxiliary
        integer(I8), allocatable    :: labels(:, :)

        ! Construct
        !==========
        allocate(labels(size(x, 1), 1))
        labels = 0
        call polygon%Construct(x, y, labels, forcevert_opt)

    end subroutine

    ! Vertex coordinate updating
    subroutine UpdatePolygonVertexCoordinates(polygon, x, y)

        ! Description
        !============
        ! This routine updates the polygon coordinates using the new
        ! values set by x, y. These values should have the same size as
        ! the original coordinates (otherwise, construct a new polygon).
        ! Metrics are recomputed as required, but duplicate points are 
        ! NOT removed anymore. If this is desired, simply reconstruct the
        ! polygon from scratch as this will be almost the same cost. 

        ! Notes:
        !=======
        ! Note 1: the advantage of not removing any vertices is that the
        ! topology of the polygon remains the same, which may be 
        ! required for shape optimization purposes, for example. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)               :: polygon 
        real(R8), intent(in)            :: x(:), y(:)

        ! Initialize
        !===========
        ! Check sizes 
        if ( (size(x, 1) /= size(y, 1)) .or. (size(x, 1) /= polygon%nv)) then 
            ! Throw error
            call gdErrorHandler('UpdatePolygonVertexCoordinates: ' // &
                'incompatible sizes of new coordinates, check input')
        end if

        ! Attribute
        !==========
        polygon%x = x 
        polygon%y = y

        ! Update fields
        !==============
        ! Check if the polygon is closed
        call polygon%IsClosedPolygon()

        ! Construct the polygon metrics
        call polygon%ComputeMetrics()

        ! Check if the polygon is simple
        call polygon%IsSimplePolygon()

        ! Check if the polygon self intersects
        ! call polygon%IsSelfIntersectingPolygon()

    end subroutine

    ! Polygon refiner
    subroutine RefinePolygonUniform(polygon, dlmax, minsplit)

        ! Description
        !============
        ! This routine refines the edges of an existing polygon by 
        ! inserting additional points at edges that are too long. 
        ! The maximal length that an edge should be is determined by 
        ! dlmax. Refinement is done in a uniform way.

        ! Additionally, a minimal refinement level can be set using 
        ! minsplit, which indicates the minimal amount of splittings 
        ! per edge that should occur. If set to zero, this has no 
        ! influence. This setting can help to split very small edges
        ! without blowing up the polygon size, or to uniformly refine
        ! the entire polygon but keeping the original resolution 
        ! distribution (then set dlmax to a very large value to have 
        ! no effect)

        ! Note: the labels of original vertices are retained, but 
        ! the new vertices' labels are simply initialized to zero since
        ! there is not an easy way to propagate these. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                       :: polygon 
        real(R8), intent(in)                    :: dlmax
        integer(I8), intent(in)                 :: minsplit

        ! Auxiliary
        integer(I8), allocatable, dimension(:,: )   :: newlabels
        integer(I8), allocatable, dimension(:)      :: nne

        real(R8), allocatable, dimension(:)     :: le, newx, newy, &
            dxe, dye, frac, dx, dy, x, y

        ! Loop
        integer(I8)                             :: i, vc, k

        ! Initialize
        !===========
        ! Unpack
        associate(&
            edges       => polygon%edges,   &
            ne          => polygon%ne,      &
            nv          => polygon%nv,      &
            labels      => polygon%labels,  &
            vert        => polygon%vert,    &
            xv          => polygon%x,       &
            yv          => polygon%y        &
            )

        ! Precompute
        x = xv(vert)
        y = yv(vert)
        dxe = x(2:ne+1) - x(1:ne)
        dye = y(2:ne+1) - y(1:ne)
        le = sqrt( dxe**2 + dye**2 )
        if (dlmax > 0.0_R8) then 
            nne = floor(le/dlmax)
        else
            nne = floor(0.0_R8*le)
        end if 
        where (nne < minsplit) nne = minsplit

        ! Initialize new edges
        allocate(newx(sum(nne)+ne+1), newy(sum(nne)+ne+1), &
            newlabels(sum(nne)+ne+1, size(labels, 2)))

        ! Construct
        !==========
        ! Initialize vertex
        vc = 0
        
        ! First coordinate
        vc = vc + 1
        newx(vc) = x(1) 
        newy(vc) = y(1)
        newlabels(vc, :) = labels(1, :)

        ! Loop over all old edges
        do i = 1, ne 

            if (nne(i) > 0) then 
                ! Compute segment coordinates
                frac = real([(k, k = 1, nne(i))], kind=R8)/real(nne(i)+1, kind=R8)
                dx = frac*dxe(i)
                dy = frac*dye(i)

                ! Add
                newx(vc+1:vc+nne(i)) = newx(vc) + dx 
                newy(vc+1:vc+nne(i)) = newy(vc) + dy 
                newlabels(vc+1:vc+nne(i), :) = 0_I8

                ! Update counter
                vc = vc + nne(i)

            end if 

            ! Last coordinate
            vc = vc + 1
            newx(vc) = x(i+1) 
            newy(vc) = y(i+1)
            newlabels(vc, :) = labels(vert(i+1), :)

        end do

        ! Construct polygon
        !==================
        call polygon%Construct(newx, newy, newlabels)

        ! Housekeeping
        end associate
        

    end subroutine

    ! Polygon compressor
    subroutine CompressPolygon(polygon)

        ! Description
        !============
        ! Wrapper of the CompressSimplePolygon for polygon objects. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)               :: polygon

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: x, y
        integer(I8), allocatable, dimension(:, :)   :: labels

        ! Loop

        ! Initialize
        !===========
        ! Unpack original vertices ()
        x = polygon%x(polygon%vert)
        y = polygon%y(polygon%vert)
        labels = polygon%labels(polygon%vert, :)

        ! Compress
        !=========
        call CompressSimplePolygon(x, y, labels)

        ! Reconstruct
        !============
        call Polygon%Construct(x, y, labels)

    end subroutine

    ! Simple polygon compression
    subroutine CompressSimplePolygon(x, y, labels)

        ! Description
        !============
        ! Compress a polygon by removing vertices that lie on straight
        ! line segments (measured by computing the angle between 
        ! neighbouring edges and comparing to pi). Machine precision 
        ! effects may occur though, so use wisely. Useful for polygons
        ! that have previously been refined and have to be coarsened 
        ! again to their initial state. Should (apart from aforementioned
        ! precision effects) not influence the shape of the polygon. 
        ! Differences in labels are not accounted for. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), allocatable, intent(inout)        :: x(:), y(:) 
        integer(I8), allocatable, intent(inout)     :: labels(:, :)

        ! Auxiliary
        integer(I8)                                 :: np
        real(R8), allocatable, dimension(:)         :: dx, dy, theta, &
            tempx, tempy
        logical, allocatable, dimension(:)          :: keepind 
        integer(I8), allocatable                    :: templabels(:, :)

        ! Loop
        integer(I8)                                 :: i 

        ! Initialize
        !===========
        np = size(x)
        allocate(keepind(np-2)) ! first and last point can never be deleted
        keepind = .true.

        ! Compute angles
        !===============
        dx = x(2:) - x(1:np-1)
        dy = y(2:) - y(1:np-1)
        theta = atan2(dx(1:np-2)*dy(2:np-1) - dx(2:np-1)*dy(1:np-2), &
            dx(1:np-2)*dx(2:np-1) + dy(2:np-1)*dy(1:np-2))

        ! Simplify
        !=========
        where (abs(theta) <= 1e-8_R8) keepind = .false.
        allocate(tempx(count(keepind)+2), tempy(count(keepind)+2), &
            templabels(count(keepind)+2, size(labels, 2)))
        tempx = [x(1), pack(x(2:size(x)-1), keepind), x(size(x))]
        tempy = [y(1), pack(y(2:size(y)-1), keepind), y(size(y))]
        do i = 1, size(labels, 2)
            templabels(:, i) = [labels(1, i), pack(labels(2:size(x)-1, i), keepind), labels(size(x), i)]
        end do 
        x = tempx
        y = tempy
        labels = templabels  

    end subroutine
    
    ! Allocator
    subroutine AllocatePolygon(polygon)

        ! Description
        !============
        ! Allocate the polygon fields for a given number of vertices nv
        ! and edges ne. These numbers must be attributed beforehand. 

        ! Declare variables
        !=================
        ! Arguments
        class(PolygonUDT)        :: polygon


        ! Allocate
        !=========
        associate(nv => polygon%nv, ne => polygon%ne, nl => polygon%nl)

        ! Check
        if (allocated(polygon%x)) then 
            ! Deallocate
            call polygon%Deallocate()
        end if 

        allocate(polygon%x(nv), polygon%y(nv), polygon%vert(ne+1), &
            polygon%edges(ne, 2), polygon%ex(ne), polygon%ey(ne), &
            polygon%tx(ne), polygon%ty(ne), polygon%tn(ne), &
            polygon%nx(ne), polygon%ny(ne), polygon%nn(ne), &
            polygon%labels(nv, nl))

        end associate


    end subroutine

    ! Deallocator
    subroutine DeallocatePolygon(polygon)

        ! Description
        !============
        ! Deallocate polygon fields

        ! Declare variables
        !=================
        ! Arguments
        class(PolygonUDT)        :: polygon

        ! Deallocate
        !===========
        deallocate(polygon%x, polygon%y, polygon%vert, &
            polygon%edges, polygon%ex, polygon%ey, &
            polygon%tx, polygon%ty, polygon%tn, &
            polygon%nx, polygon%ny, polygon%nn, &
            polygon%labels)


    end subroutine

    ! Initialization
    subroutine InitializePolygon(polygon, x, y, labels)

        ! Description
        !============
        ! Initialize the polygon by allocating the required fields, 
        ! computing the initial polygon dimensions, and adding the 
        ! coordinates. It is assumed that all quantities are compliant
        ! with the polygon definition. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 
        real(R8), intent(in)                :: x(:), y(:)
        integer(I8), intent(in)             :: labels(:, :)

        ! Auxiliary

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Sizes
        polygon%nv = size(x)
        polygon%ne = polygon%nv - 1
        polygon%nl = size(labels, 2)

        ! Allocate
        call polygon%Allocate()

        ! Set coordinates
        polygon%x = x 
        polygon%y = y 

        ! Set edges
        polygon%edges(:, 1) = [(k, k = 1, polygon%nv-1)]
        polygon%edges(:, 2) = [(k, k = 2, polygon%nv)]

        ! Set labels (currently following x, y)
        polygon%labels = labels


    end subroutine

    ! Duplicate point remover
    subroutine RemoveDuplicatePoints(polygon, forcevert_opt)

        ! Description
        !============
        ! This function checks whether coordinates of a polygon appear multiple
        ! times and replaces them by a single coordinate (and renumbers the
        ! vertices accordingly). If the polygon is destroyed in the process (i.e.
        ! all points are duplicate), an error is thrown.
        ! Note that this is an expensive routine
        ! and should only be called once during construction of the polygon! 

        ! Input
        !------
        ! - polygon:        polygon structure with the x, y, vert, labels fields
        ! - forcevert_opt:  (optional) force vertices to be in polygon 
        !                   (i.e. don't allow them to be deleted). This 
        !                   may result in unexpected outcomes, warnings 
        !                   will be printed in that case

        ! Output
        !-------
        ! - polygon:        same structure, but with duplicate points removed and
        !                   empty polygons destroyed.

        ! Algorithm
        !----------
        ! Loop over all points, compute the distance w.r.t. each consecutive point
        ! in the polygon. If the distance is smaller than a small tolerance, the
        ! point is considered to be duplicate. In that case, the next point is
        ! adapted to be the current point, the vertices are changed, and the
        ! algorithm continues on the next points. 

        ! Notes 
        !------
        ! Note 1: the polygon metrics are *not* updated in this routine - this has
        ! to be done afterwards!

        ! Note 2: when duplicate points are removed, the polygon has to 
        ! be rebuilt and dimensions of arrays will change. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)           :: polygon 
        logical, intent(in), optional   :: forcevert_opt(:)

        ! Auxiliary
        real(R8)                    :: d 

        real(R8), allocatable       :: tempx(:), tempy(:)

        integer(I8)                 :: nlabels
        integer(I8), allocatable    :: mapping(:), tempedge(:, :), &
            diffindex(:), templabels(:, :)
        logical                     :: vertnotremoved
        logical, allocatable        :: keepvertex(:), keepedges(:), &
            forcevert(:)

        ! Loop
        integer(I8)                 :: i, j, k

        ! Find duplicate points
        !======================
        ! Associate
        associate(&
            x       => polygon%x,       &
            y       => polygon%y,       &
            nv      => polygon%nv,      &
            ne      => polygon%ne,      & 
            edges   => polygon%edges)

        ! Initialize
        allocate(mapping(nv), keepvertex(nv), diffindex(nv)) 
        mapping(:) = [(k, k = 1, nv)]
        keepvertex(:) = .true.
        diffindex(:) = 0
        nlabels = size(polygon%labels, 2)
        if (present(forcevert_opt)) then 
            ! Assign
            forcevert = forcevert_opt

            ! Checks
            if (size(forcevert) /= nv) then 
                call PolygonErrorHandler('RemoveDuplicatePoints: the ' // & 
                    'forcevert optional input argument should have the ' // & 
                    'same dimensions as the polygon coordinates, check input')
            end if
        else
            allocate(forcevert(nv))
            forcevert = .false. 
        end if 
        vertnotremoved = .false. 
        
        ! Loop over all points but one 
        do i = 1, nv-1 
            ! Compute distance to the next coordinates
            do j = i+1, nv 
                ! Compute the distance
                call Distance(d, x(i), y(i), x(j), y(j))

                ! Check if distance is smaller than distance tolerance
                if (d < disttol) then 
                    ! Mark for removal
                    if (.not. forcevert(j)) then 
                        ! Take vertex j
                        keepvertex(j) = .false. 

                        ! Update mapping
                        mapping(j) = mapping(i) ! not just i - i can already be mapped to another vertex
                    else
                        ! Check if we can remove vertex i
                        if (.not. forcevert(i)) then 
                            ! Take vertex i 
                            keepvertex(i) = .false.

                            ! Update mapping
                            mapping(i) = mapping(j)
                        else
                            ! No vertex can be removed, but issue warning
                            vertnotremoved = .true. 
                        end if 
                    end if 
                end if
            end do
        end do

        ! Stop associate
        end associate

        ! Issue warning if necessary
        if (vertnotremoved) then 
            call PolygonWarningHandler('RemoveDuplicatePoints: vertices ' &
                // 'were duplicate but not removed by forcing. Results may be ' & 
                // 'surprising')
        end if 

        ! Remove & rebuild
        !=================
        ! Compute offset for new vertices (renumbering due to deletion)
        where (.not. keepvertex) diffindex = 1
        do i = 2, polygon%nv 
            diffindex(i) = diffindex(i-1) + diffindex(i)
        end do

        ! Remap the vertices in the edges
        polygon%edges(:, 1) = mapping(polygon%edges(:, 1))
        polygon%edges(:, 2) = mapping(polygon%edges(:, 2))
        polygon%edges(:, 1) = polygon%edges(:, 1) - diffindex(polygon%edges(:, 1))
        polygon%edges(:, 2) = polygon%edges(:, 2) - diffindex(polygon%edges(:, 2))


        ! Check if any edges have to be removed
        keepedges = polygon%edges(:, 1) .ne. polygon%edges(:, 2)
        if (.not. all(keepedges)) then 
            call PolygonWarningHandler('RemoveDuplicatePoints: edges ' &
                // 'removed that consist of only duplicate points')
        end if 

        ! Check if the polygon is destroyed in the process
        if (count(keepvertex) < 2) then 
            call PolygonErrorHandler('RemoveDuplicatePoints: removal ' &
                // 'of duplicate points resulted in polygon ' &
                // 'destruction - check input.')
        end if 

        ! Rebuild the polygon
        allocate(tempedge(count(keepedges), 2), &
            tempx(count(keepvertex)), tempy(count(keepvertex)), &
            templabels(count(keepvertex), nlabels))
        tempedge(:, 1) = pack(polygon%edges(:, 1), keepedges)
        tempedge(:, 2) = pack(polygon%edges(:, 2), keepedges)
        tempx = pack(polygon%x, keepvertex)
        tempy = pack(polygon%y, keepvertex)
        do i = 1, nlabels
            templabels(:, i) = pack(polygon%labels(:, i), keepvertex)
        end do

        call polygon%Deallocate() 
        polygon%nv = count(keepvertex)
        polygon%ne = count(keepedges)
        call polygon%Allocate() 
        polygon%edges = tempedge 
        polygon%x = tempx 
        polygon%y = tempy
        polygon%labels = templabels

        ! Housekeeping
        !=============
        deallocate(tempedge, tempx, tempy, diffindex, keepvertex, &
            keepedges, templabels)


    end subroutine

    ! Vertex setter
    subroutine SetPolygonVertices(polygon)

        ! Description
        !============
        ! This routine extracts the subsequent set of vertices of a 
        ! polygon from its edges. This set of vertices is then stored
        ! in the polygon%vert field. Note that the size of this field
        ! is not necessarily equal to the size of polygon%x! This is 
        ! because multiple points can occur multiple times in the
        ! polygon. Note that labels still follow the x, y sequence

        ! It is assumed that the polygon has been properly initialized
        ! and that any duplicate points have been removed (will work 
        ! just as well without this last requirement, but then vert will
        ! change after removing duplicate points). 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 

        ! Auxiliary
        integer(I8)                         :: cvind, nv  
        logical                             :: check(1:4)

        ! Loop
        integer(I8)                         :: k 

        ! Initialize
        !===========
        ! Associate
        associate( & 
            ne          => polygon%ne,      &
            edges       => polygon%edges,   &
            vert        => polygon%vert) 
        

        ! Check
        check(:) = .false. 

        ! Construct vert
        !===============
        ! Hedge for single edge polygon
        if (ne == 1) then 
            ! Sorted by default, return
            vert(1) = edges(1, 1)
            vert(2) = edges(1, 2)
            return 
        end if 

        ! Get the starting edge
        check(1) = edges(1, 1) == edges(2, 1)
        check(2) = edges(1, 1) == edges(2, 2)
        check(3) = edges(1, 2) == edges(2, 1)
        check(4) = edges(1, 2) == edges(2, 2)

        ! Do checks
        if ( (check(1) .and. ( .not. any(check(2:4))) ) .or. &
                (check(2) .and. ( .not. any(check([1, 3, 4]))) ) ) then 
            ! Second vertex is starting vertex
            cvind = 2
        elseif ( (check(4) .and. ( .not. any(check(1:3))) ) .or. &
                (check(3) .and. ( .not. any(check([1, 2, 4]))) ) ) then 
            ! First vertex is starting vertex
            cvind = 1
        else
            ! Something wrong with input - throw error
            call PolygonErrorHandler('ExtractPolygonVertices: ' &
                // 'something wrong with input - please check input ' &
                // 'polygon')
        end if

        ! Loop 
        k = 1
        do while (k < ne) 
            ! Get the current vertex and add
            vert(k) = edges(k, cvind) 

            ! Get the next vertex
            if (cvind == 1) then 
                nv = edges(k, 2)
            else 
                nv = edges(k, 1)
            end if 

            ! Update cvind 
            if (nv == edges(k+1, 1)) then 
                cvind = 1
            elseif (nv == edges(k+1, 2)) then 
                cvind = 2
            else 
                ! Throw error, next edge does not contain the vertex
                call PolygonErrorHandler('ExtractPolygonVertices: ' &
                    // 'could not find next edge - check input polygon')
            end if 

            ! Update counter
            k = k + 1
        end do 
 
        ! Add last vertices
        vert(k) = edges(ne, cvind)
        if (cvind == 1) then 
            vert(k+1) = edges(ne, 2)
        else 
            vert(k+1) = edges(ne, 1)
        end if 

        end associate 

    end subroutine

    ! Vertex ID getter
    subroutine GetPolygonVertexID(polygon, ID)

        ! Description
        !============
        ! Get polygon vertex ID vector - this is simply the 1:nv 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                       :: polygon 
        integer(I8), allocatable, intent(out)   :: ID(:)

        ! Loop
        integer(I8)                             :: k 

        ! Add
        !====
        ID = [(k, k = 1, polygon%nv)]

    end subroutine

    ! Flipper
    subroutine FlipPolygon(polygon)

        ! Description
        !============
        ! This routine 'flips', i.e. changes the direction, of a polygon 
        ! by changing the order of the vertices (not the vertex 
        ! coordinates!). The metrics do not change, except for the sign 
        ! of the normal and tangent vectors. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)       :: polygon 

        ! Flip
        !=====
        ! Associate
        associate( &
            nv    => polygon%nv, &
            ne    => polygon%ne)

        ! Vertices
        polygon%vert    = polygon%vert(ne+1:1:-1)
        polygon%edges   = polygon%edges(ne:1:-1, 2:1:-1)
        
        ! Edge centers
        polygon%ex      = polygon%ex(ne:1:-1)
        polygon%ey      = polygon%ey(ne:1:-1)

        ! Tangents and normals
        polygon%tx      = -polygon%tx(ne:1:-1)
        polygon%ty      = -polygon%ty(ne:1:-1)
        polygon%tn      = polygon%tn(ne:1:-1)
        polygon%nx      = -polygon%nx(ne:1:-1)
        polygon%ny      = -polygon%ny(ne:1:-1)
        polygon%nn      = polygon%nn(ne:1:-1)

        ! Housekeeping
        !=============
        end associate 
        

    end subroutine

    ! Is polygon closed
    subroutine IsClosedPolygon(polygon)

        ! Description
        !============
        ! Check if the polygon is closed upon itself. It is assumed that the
        ! polygon's duplicate vertices are removed and that as such simply the
        ! start and end edge can be checked for a duplicate vertex. If one (or more
        ! in the case of a polygon with two edges with different orientation) are
        ! found, the polygon is said to be closed. This is then updated in the
        ! polygon%isclosed logical field.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)               :: polygon

        ! Checks
        !=======
        ! Simply check start and end vertex
        if (polygon%vert(1) == polygon%vert(polygon%ne+1)) then 
            polygon%isclosed = .true. 
        else 
            polygon%isclosed = .false.
        end if 

    end subroutine

    ! Metrics computation
    subroutine ComputePolygonMetrics(polygon)

        ! Description
        !============
        ! Compute the geometric quantities (or metrics) of the polygon. It is
        ! assumed that the polygon's duplicate points have been removed (otherwise
        ! it will work just as well, but NaNs may occur in the tangents/normals). 

        ! Algorithm
        !==========
        ! Most things are quite straightforward computations. To determine the
        ! inward pointing normal, we construct for each edge pair a point projected
        ! along the current normal and determine whether point lies within the 
        ! polygonal region using 'inpolygon'. If this is not the case, the sign of
        ! the normal is reversed. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                   :: polygon 

        ! Auxiliary
        real(R8)                            :: testdist
        real(R8), allocatable, dimension(:) :: xe1, xe2, ye1, ye2, &
            temptn, tpx, tpy

        logical, allocatable                :: isin(:)

        ! Loop
        integer(I8)                         :: i

        ! Initialize
        !===========
        ! Distance to construct test points for normal
        testdist = disttol*100 ! larger than tolerance 

        ! Check to be sure if the polygon is closed
        call polygon%IsClosedPolygon()

        ! Compute metrics
        !================
        ! Associate
        associate(&
            x           => polygon%x,       &
            y           => polygon%y,       &
            ne          => polygon%ne,      &
            tx          => polygon%tx,      &
            ty          => polygon%ty,      &
            tn          => polygon%tn,      &
            nx          => polygon%nx,      &
            ny          => polygon%ny,      &
            nn          => polygon%nn,      &
            ex          => polygon%ex,      &
            ey          => polygon%ey,      &
            edges       => polygon%edges,   &
            vert        => polygon%vert)

        ! Polygon edge coordinates
        allocate(xe1(ne), ye1(ne), xe2(ne), ye2(ne), temptn(ne))
        xe1 = x(edges(:, 1)) 
        xe2 = x(edges(:, 2))
        ye1 = y(edges(:, 1)) 
        ye2 = y(edges(:, 2))

        ! Edge centers
        ex = 0.5*(xe1 + xe2)
        ey = 0.5*(ye1 + ye2)

        ! Tangents
        tx = xe2 - xe1
        ty = ye2- ye1
        call Distance(temptn, xe1, ye1, xe2, ye2)
        tn = temptn

        ! Initial normals
        nx = -ty
        ny = tx
        nn = tn 

        ! Determine normal orientation
        !=============================
        ! Construct test points
        allocate(tpx(ne), tpy(ne))
        tpx = ex + testdist*nx/nn 
        tpy = ey + testdist*ny/nn 

        ! See if points lie in polygon
        call polygon%Inpolygon(tpx, tpy, isin)
        where (.not. isin) 
            ! Switch sign
            nx = -nx 
            ny = -ny
            tx = -tx 
            ty = -ty
        end where

        ! Don't forget to switch edges!
        do i = 1, ne
            if (.not.isin(i)) then 
                polygon%edges(i, :) = polygon%edges(i, [2, 1])
            end if 
        end do


        end associate

        ! Housekeeping
        !=============
        deallocate(xe1, ye1, xe2, ye2, tpx, tpy, temptn)

    end subroutine

    ! Simple polygon determination
    subroutine IsSimplePolygon(polygon)

        ! Description
        !============
        ! Check if the polygon is a simple polygon, meaning that each vertex only
        ! belongs to maximally two edges. It is assumed that duplicate vertices and
        ! pathological edges have already been removed, and that therefore this
        ! check can be done on a vertex ID basis only. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)           :: polygon 

        ! Auxiliary
        integer(I8), allocatable    :: counter(:) 

        ! Loop
        integer(I8)                 :: i


        ! Initialize
        !===========
        allocate(counter(size(polygon%x))) 

        ! Check
        !======
        ! Count vertex occurrences 
        counter(:) = 0
        do i = 1, size(polygon%edges, 1) 
            counter(polygon%edges(i, 1)) = counter(polygon%edges(i, 1)) + 1 
            counter(polygon%edges(i, 2)) = counter(polygon%edges(i, 2)) + 1 
        end do 

        ! Set logical
        polygon%issimple = all(counter <= 2)

        ! Housekeeping
        !=============
        deallocate(counter) 



    end subroutine

    ! Self-intersections of a polygon
    function IsSelfIntersectingPolygon(polygon) result(isselfintersecting)

        ! Description
        !============
        ! This routine sets the 'isselfintersecting' field by computing
        ! whether any self intersections are present. Non-simple 
        ! polygons are by definition self-intersecting, closed 
        ! polygons are not. It is assumed all necessary metrics have been 
        ! computed and that all other fields (issimple, isclosed) are 
        ! up to date. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)               :: polygon 
        logical                         :: isselfintersecting

        ! Auxiliary
        real(R8), allocatable           :: x(:), y(:) 
        integer(I8), allocatable        :: s1(:), s2(:)

        ! Loop

        ! Initialize
        !===========
        isselfintersecting = .false. 

        ! Checks
        !=======
        ! Need to check self intersections
        call polygon%SelfIntersections(x, y, s1, s2)

        if (size(x) > 0) then 
            ! Self-intersections found, need to check 
            if (size(x) > 1) then 
                isselfintersecting = .true. 
            else 
                ! Additional check on closure
                if (polygon%isclosed) then 
                    if ( (s1(1) .ne. 1) .and. (s2(1) .ne. 1) ) then
                        ! Something wrong here, this shouldn't be happening
                        call PolygonErrorHandler('SelfIntersections: ' &
                            // 'closed polygon with single ' &
                            // 'intersection that is not on end ' &
                            // 'points detected - chck input')
                    end if 
                    ! Otherwise, do nothing - closed polygons are 
                    ! not considered self-intersecting
                else
                    isselfintersecting = .true. 
                end if 
            end if
        end if 

    end function

    !------------------------------------------------------------------!
    !                     Polygon derivative routines                  !
    !------------------------------------------------------------------!
        

    !------------------------------------------------------------------!
    !                              Auxiliary                           !
    !------------------------------------------------------------------!

    ! Error handler
    subroutine PolygonErrorHandler(msg) 

        ! Description
        !============
        ! Error handler. Will print the message in 'msg' and stop the
        ! execution of the program. 

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: msg 

        ! Handle error
        !=============
        ! Print 
        print *, 'PolygonErrorHandler: ', msg 

        ! Stop program execution
        stop 

    end subroutine

    ! Warning handler
    subroutine PolygonWarningHandler(msg) 

        ! Description
        !============
        ! Warning handler: will print only the message with warning in
        ! front, will not stop program execution.

        ! Declare variables
        !==================
        ! Arguments
        character(*), intent(in)        :: msg 

        ! Handle error
        !=============
        ! Print 
        print *, 'Warning: PolygonWarningHandler: ', msg 

    end subroutine

    ! Polygon edge sorter
    subroutine SortPolygonEdges(pein, ne, sortindex, ispolygonstart, &
            isbranchingpolygon, polygonID)

        ! Description
        !============
        ! This routine sorts the edges of a simply polygon. The definition
        ! that is used here for simple polygon demands that the polygon 
        ! either closes perfectly on itself, or is a non-branching polygon 
        ! (i.e. each vertex only has maximal two edges). This is checked 
        ! during the routine execution when determining the polygon edges. 
        ! Multiple open and closed polygons are supported. The logical 
        ! 'ispolygonstart' indicates which of the (sorted!) edges is the 
        ! start of a new polygon. 

        ! Update: added support for branching polygons. These are 
        ! returned as separate polygon parts. If one wants to check if
        ! it is a branching polygon or not, one can check if there are any true
        ! indices in 'isbranchingpolygon'. Note that the branching polygons are
        ! represented by an ensemble of non-branching polygons and can therefore be
        ! reconstructed if necessary. Retrieving all parts of a branching polygon
        ! has to be done 'manually' by checking for each branching polygon which
        ! vertices it has in common with another one. Support for this might be
        ! added in the future. 
        ! Note: the current implementation to do this should have 
        ! the same basic mechanics as finding the subgraphs of an 
        ! undirected graph. Perhaps in the future, a graph-based 
        ! implementation could be directly used here to increase 
        ! readability and maintainability... 
    
        ! Arguments
        !==========
        !
        ! - pein:           (input) ne-by-2 array of vertex indices that has
        !                   to be sorted. 
        ! - ne:             (input) number of edges (integer)
        ! - sortindex:      (output) index that contains the sequence of 
        !                   sorted polygon edge indices, i.e. 
        !                   pein(sortindex,:) should sort the edges.
        ! - ispolygonstart: (output) ne-by-1 logical of which the i'th
        !                   element is true if sortindex(i) is the start 
        !                   index of a new polygon. The edges of this 
        !                   polygon are all edges between this true value 
        !                   and the next. 
        ! - isbranchingpolygon  : ne-by-1 array indicating if the polygon
        !                       is branching or not
        ! - polygonID           :: ne-by-1 array with polygon IDs 
        !                       (optional). Useful to reconstruct branching
        !                       polygons
    
        ! Algorithm
        !==========
        !
        ! 0)    Initialize and allocate
        ! 1)    Find a starting vertex:
        !
        !       Take a vertex, check how many times it occurs.
        ! 
        !       If one: 
        !       start vertex found, go to 2). Else, repeat 1). If no
        !       vertex found (and has not errored meanwhile), only closed
        !       polygons are left.
        !       Take any vertex as starting vertex, go to 2)
        !
        !       If two:
        !       inner vertex found of a polygon, repeat 1) for the next
        !       vertex.
        !
        !       If more than two:
        !       branching polygon, this vertex has to be a starting 
        !       (or ending) vertex. 
        !        
        ! 2)    Find the edges of the current polygon. 
        !
        !       2.1)    Find the edge that contains the current vertex which 
        !               initially is the starting vertex and which has not
        !               yet been sorted. If no edges meet this criterion,
        !               exit the loop and go to 3). If one edge is found,
        !               go to 2.2). If multiple edges are found, throw error
        !               (this should be impossible though and would indicate
        !               a bug in the code)
        !       2.2)    For the next edge, find the other vertex and set it 
        !               to the current vertex. Go to 2.1)
    
        ! 3)    Update data structures, and check if any edges remain. If 
        !       yes, go to 1), if no, exit the routine and return. 
    
        ! Initialize
        !===========
        ! The usual
        implicit none
    
        ! Declare variables
        !==================
        ! Input
        integer(I8), dimension(ne,1:2), intent(in)  :: pein ! polygon edges 
        integer(I8), intent(in)                     :: ne
        
        ! Output
        integer(I8), dimension(ne)      :: sortindex 
        logical, dimension(ne)          :: ispolygonstart, isbranchingpolygon
        integer(I8), allocatable, optional  :: polygonID(:)
    
        ! Mixed
    
        ! Loop
        logical                 :: allfound, startfound, polygonfound
        integer(I8)             :: i, j, k, spind
    
        ! Auxiliary
        integer(I8)                 :: nv, nremedges, nextremedge
    
        logical                     :: allbranchingfound
        logical, allocatable        :: isedgesorted(:), isremedgesorted(:), &
            mask(:), isbranchingvertex(:, :), remisbranching(:, :), &
            hasbv(:), sortedisbranchingvertex(:, :), &
            isnonbranchingstartvertex(:, :), isvertfound(:), tracevert(:)
    
        integer(I8)                 :: pID, ind 
        integer(I8), allocatable    :: remedges(:,:), edgeID(:), &
            remedgeID(:), temparray(:), allv(:), oc(:), el(:), sortind(:), &
            alloc(:), ps(:), tv(:), tvu(:), tvID(:), sortededges(:, :), &
            pe(:), pse(:, :), pee(:, :), sortind2(:)
    
        ! Main program
        !=============
        ! Check
        if (size(pein,2) /= 2) then
            ! Throw error
            call gdErrorHandler('SortPolygonEdges: input argument pein should be a ne-by-2 integer array')
    
        end if

        if (size(pein, 1) == 0) then 
            ! Simply return, already sorted
            return
        end if
    
        ! Initialize
        pID = 1
        allocate(isedgesorted(ne)) ! logical to indicate if edge has been sorted
        allocate(edgeID(ne))
    
        ispolygonstart = .false.
        isedgesorted = .false.
        isbranchingpolygon = .false. 
        edgeID(:) = (/ (i, i=1,ne,1) /)
        allfound = .false. ! while loop variable
        spind = 1 ! sorted polygon index
        sortindex = 0

        ! Check if branching polygons exist by counting occurrence
        allv = [pein(:, 1), pein(:, 2)]
        call CountOccurrence(allv, oc, el, sortindoc=sortind, sortindel=sortind2)
        alloc = oc(sortind)
        allocate(isbranchingvertex(ne, 2), isnonbranchingstartvertex(ne, 2))
        isbranchingvertex(:, 1) = alloc(1:ne) > 2_I8
        isbranchingvertex(:, 2) = alloc(ne+1:2*ne) > 2_I8 
        isnonbranchingstartvertex(:, 1) = alloc(1:ne) == 1_I8
        isnonbranchingstartvertex(:, 2) = alloc(ne+1:2*ne) == 1_I8
    
        ! Loop
        allbranchingfound = count(isbranchingvertex) == 0
        do while (allfound .eqv. .false.)
            ! Set the polygon starting index
            ispolygonstart(spind) = .true.
    
            ! Allocate inner loop variables
            nremedges = count(isedgesorted .eqv. .false.)
            allocate(remedges(nremedges,2))
            allocate(remedgeID(nremedges))
            allocate(isremedgesorted(nremedges))
            allocate(mask(nremedges))
            allocate(remisbranching(nremedges, 2))
    
            ! Get the remaining edges
            remedges(:,1) = pack(pein(:,1), (isedgesorted .eqv. .false.))
            remedges(:,2) = pack(pein(:,2), (isedgesorted .eqv. .false.))
            remedgeID(:) = pack(edgeID, (isedgesorted .eqv. .false.))
            remisbranching(:, 1) = pack(isbranchingvertex(:, 1), .not. isedgesorted)
            remisbranching(:, 2) = pack(isbranchingvertex(:, 2), .not. isedgesorted)
            isremedgesorted(:) = .false.
    
            ! Find a starting vertex of a branching vertex
            startfound = .false. 
            k = 1
            do while ((startfound .eqv. .false.) .and. (k <= nremedges) .and. &
                .not. allbranchingfound)

                ! Check if there are any branching edges among the 
                ! remaining edges
                if (remisbranching(k, 1)) then 

                    ! Found starting point
                    startfound = .true. 
                    nv = remedges(k, 2)

                    ! Set as branching polygon
                    isbranchingpolygon(remedgeID(k)) = .true.
                    
                    ! Add the current edge
                    sortindex(spind) = remedgeID(k)
    
                    ! Set edge as sorted
                    isremedgesorted(k) = .true.
    
                    ! Update indices
                    k = k+1
                    spind = spind+1

                elseif (remisbranching(k, 2)) then 

                    ! Found starting point
                    startfound = .true. 
                    nv = remedges(k, 1)

                    ! Set as branching polygon
                    isbranchingpolygon(remedgeID(k)) = .true.
                    
                    ! Add the current edge
                    sortindex(spind) = remedgeID(k)
    
                    ! Set edge as sorted
                    isremedgesorted(k) = .true.
    
                    ! Update indices
                    k = k+1
                    spind = spind+1

                else 
                    ! Next edge
                    k = k + 1
                end if

                ! Check
                if (.not. startfound .and. k >= nremedges) then 
                    ! Found all branching polygons
                    allbranchingfound = .true. 
                end if 
            end do 

            ! If no starting vertex was found, find the starting vertex 
            ! of a non-branching polygon
            k = 1
            do while ((startfound .eqv. .false.) .and. (k <= nremedges))

                ! Count how many times the current edge vertices occur
                !tc1 = count(remedges(:,1) == remedges(k,1)) & 
                !    + count(remedges(:,2) == remedges(k,1))
                !tc2 = count(remedges(:,1) == remedges(k,2)) & 
                !    + count(remedges(:,2) == remedges(k,2))
    
                ! Check & add
                if (isnonbranchingstartvertex(remedgeID(k), 1)) then
                    ! Found a start vertex
                    startfound = .true.
                    nv = remedges(k,2) ! next vertex
    
                    ! Add the current edge
                    sortindex(spind) = remedgeID(k)
    
                    ! Set edge as sorted
                    isremedgesorted(k) = .true.
    
                    ! Update indices
                    k = k+1
                    spind = spind+1
                else if (isnonbranchingstartvertex(remedgeID(k), 2)) then
                    ! Found a start vertex
                    startfound = .true.
                    nv = remedges(k,1) ! next vertex
    
                    ! Add the current edge
                    sortindex(spind) = remedgeID(k)
    
                    ! Set edge as sorted
                    isremedgesorted(k) = .true.
    
                    ! Update counters
                    k = k+1
                    spind = spind+1
                else if (isbranchingvertex(remedgeID(k), 1) .or. isbranchingvertex(remedgeID(k), 2)) then
                    ! Polygon branches, throw error
                    call gdErrorHandler('SortPolygonEdges: branching polygon detected, not supported')
                else
                    ! Next edge, update counter
                    k = k+1
                end if
            end do
    
            ! Hedge for closed polygon(s)
            if (startfound .eqv. .false.) then
                ! No start index found, yet edges remain -> has to be closed
                ! polygon. Simply take the first vertex of the first edge
                startfound = .true.
                nv = remedges(1,2)
    
                ! Add the current edge
                sortindex(spind) = remedgeID(1)
    
                ! Set edge as sorted
                isremedgesorted(1) = .true.
    
                ! Update counters
                spind = spind+1
            end if
    
            ! Find the edges of this polygon
            polygonfound = .false.
            do while (polygonfound .eqv. .false.)
                ! The next edge has the current vertex and is not yet found.
                mask(:) = (isremedgesorted .eqv. .false.) ! first requirement
                mask = mask .and. & 
                    ((remedges(:,1) == nv) .or. (remedges(:,2) == nv)) ! second requirement
                
                ! Check if this is the final edge
                if (count(mask) == 0) then
                    ! All edges were found, exit
                    polygonfound = .true. 
                else if (count(mask) > 1) then
                    ! The current vertx is a branching vertex, exit
                    polygonfound = .true.
                else
                    ! Get the next edge
                    allocate(temparray(1)) ! avoid rank conflicts
                    temparray = pack((/ (i, i=1,nremedges,1) /),mask)
                    nextremedge = temparray(1)
    
                    ! Add edge
                    sortindex(spind) = remedgeID(nextremedge)
    
                    ! Set edge as sorted
                    isremedgesorted(nextremedge) = .true.
    
                    ! Get next vertex & check if it is branching. If so, exit
                    if (remedges(nextremedge,1) == nv) then
                        nv = remedges(nextremedge,2)
                        if (remisbranching(nextremedge, 2)) then 
                            polygonfound = .true. 
                        end if 
                    else
                        nv = remedges(nextremedge,1)
                        if (remisbranching(nextremedge, 1)) then 
                            polygonfound = .true. 
                        end if 
                    end if
                    
                    ! Update counters
                    spind = spind+1

                    ! Deallocate
                    deallocate(temparray)
                end if 
            end do
    
            ! Update logicals
            isedgesorted(pack(remedgeID,isremedgesorted)) = .true. 
    
            ! Check if all edges have been found
            if (spind == ne+1) then ! one ahead due to updating rule here
                allfound = .true.
            end if
    
            ! Deallocate inner loop variables
            deallocate(remedges)
            deallocate(remedgeID)
            deallocate(isremedgesorted)
            deallocate(mask)
            deallocate(remisbranching)
    
        end do

        ! Check if we need to determine polygon IDs
        if (present(polygonID)) then 
            ! Allocate
            if (allocated(polygonID)) then 
                deallocate(polygonID)
            end if 
            allocate(polygonID(ne))
            polygonID = 0

            ! Check
            if (.not. any(isbranchingpolygon)) then 
                ! Easy
                do i = 1, ne
                    if (ispolygonstart(i)) then 
                        polygonID(i) = pID 
                        pID = pID + 1
                    end if 
                end do 
            else
                ! Not so easy - need to check polygon vertices etc. 
                ! Luckily, we know that each branching polygon part 
                ! should start in a branching polygon vertex

                ! Sort the edges 
                allocate(sortededges(ne, 2))
                sortededges(:, 1) = pein(sortindex, 1)
                sortededges(:, 2) = pein(sortindex, 2)
                sortedisbranchingvertex = isbranchingvertex
                sortedisbranchingvertex(:, 1) = isbranchingvertex(sortindex, 1)
                sortedisbranchingvertex(:, 2) = isbranchingvertex(sortindex, 2)
                
                ! Get all polygon start/end indices
                allocate(ps(count(ispolygonstart)))
                ps = pack([(k, k = 1, ne)], ispolygonstart)
                pe = [ps(2:size(ps))-1, ne]

                ! Get starting and ending edges
                allocate(pse(count(ispolygonstart), 2), &
                    pee(count(ispolygonstart), 2)) 
                pse(:, 1) = sortededges(ps, 1)
                pse(:, 2) = sortededges(ps, 2)
                pee(:, 1) = sortededges(pe, 1)
                pee(:, 2) = sortededges(pe, 2)

                ! Get the unique IDs of polygon branching vertices
                tv = [pack(pein(:, 1), isbranchingvertex(:, 1)), &
                    pack(pein(:, 2), isbranchingvertex(:, 2))]
                call Unique(tv, tvu)

                ! Initialize the branching vertex ID
                allocate(tvID(size(tvu)))
                tvID = 0 

                ! Start 'tracing' - keep looping until all start vertices were
                ! found 
                allocate(isvertfound(size(tvu)), tracevert(size(tvu)))
                isvertfound = .false. 
                tracevert = .false. 
                do while (.true.)
                    ! Take a new vertex
                    ind = findloc(isvertfound, .false., 1, back=.false.)

                    ! Check exit condition
                    if (ind == 0) then 
                        ! All should be found, exit
                        exit 
                    end if 

                    ! Sanity check
                    if (tvID(ind) /= 0) then 
                        ! Starting at a vertex that was already 
                        ! considered - this shouldn't happen
                        call gdErrorHandler('SortPolygonEdges: starting at ' // &
                            'a vertex that has already been added previously, ' // & 
                            'this is a bug')
                    end if 

                    ! Set ID
                    tvID(ind) = pID 

                    ! Set vertex to be traced
                    tracevert(ind) = .true. 

                    ! Keep tracing until no more tracing vertices are found
                    do while (.true.)
                        ! Get the next vertex to trace
                        ind = findloc(tracevert, .true., 1)
                        if (ind == 0) then 
                            exit 
                        end if 

                        ! Check which polygons have this vertex
                        hasbv = pse(:, 1) == tvu(ind) .or. pse(:, 2) == tvu(ind) &
                            .or. pee(:, 1) == tvu(ind) .or. pee(:, 2) == tvu(ind)
                        do j = 1, size(ps)
                            if (hasbv(j)) then 
                                ! Consistency check
                                if (polygonID(ps(j)) == 0 .or. polygonID(ps(j)) == pID) then 
                                    ! Add
                                    polygonID(ps(j)) = pID

                                    ! Check if any other vertices are 
                                    ! polygon start vertices and should
                                    ! be traced
                                    tracevert = tracevert .or. &
                                        (((pse(j, 1) == tvu) .or. (pse(j, 2) == tvu) .or. &
                                        (pee(j, 1) == tvu) .or. (pee(j, 2) == tvu)) &
                                        .and. .not. isvertfound)
                                else
                                    ! Something weird going wrong here
                                    call gdErrorHandler('SortPolygonEdges: ' // & 
                                        'different branching polygons seem to ' // & 
                                        'have common vertices, unexpected')
                                end if 
                            end if 
                        end do 

                        ! Mark the current vertex as traced and found
                        isvertfound(ind) = .true.
                        tracevert(ind) = .false.
                    end do 

                    ! Update branching polygon counter
                    pID = pID + 1
                end do 

                ! Set polygon IDs for non-branching polygons
                do i = 1, ne
                    if (ispolygonstart(i) .and. polygonID(i) == 0) then 
                        polygonID(i) = pID
                        pID = pID + 1
                    end if 
                end do

                
            end if 
        end if 
    
    end subroutine

    ! Vertex extractor from sorted edges
    subroutine ExtractPolygonVertices(pe, ne, pv)

        ! Description
        !===========
        ! This routine extracts the vertices of a polygon from a set of 
        ! sorted polygon edges which are given by their vertex IDs in the 
        ! 'pe' (polygon edges) ne-by-2 array. Only single polygons (can be 
        ! closed or open) are supported, i.e. no branching or multiple 
        ! polygons. This is checked for during the routine. 
    
        ! Arguments
        !==========
        !
        ! - pe:             (input) ne-by-2 array of vertex indices that has
        !                   to be sorted. 
        ! - ne:             (input) number of edges (integer)
        ! - pv:             (output) ne+1-by-1 array of polygon vertices
    
        ! Algorithm
        !==========
        ! 0)    Initialize & check
        ! 1)    Take the first edge
        !
        ! 2)    Check which vertex is the start vertex by comparing with the 
        !       vertices of the next edge. Set this as the next vertex. 
        !       If (the vertex is not found) then
        !           throw error
        !       else 
        !           add the current vertex
        !           set the other vertex as the next vertex
        !       end if
        !
        ! 3)    for the remaining edges, repeat 2) but with the start vertex
        !       equal to the next vertex. If the last vertex is reached, add 
        !       the remaining vertex. 
    
        ! Declare variables
        !==================
        ! Input
        integer(I8), dimension(:, :), intent(in)  :: pe ! polygon edges 
        integer                         :: ne
        
        ! Output
        integer(I8), allocatable, dimension(:), intent(out)    :: pv
    
        ! Mixed
    
        ! Loop
        integer(I8)                     :: k, cvind, nv
    
        ! Auxiliary
        logical, dimension(1:4)         :: check
    
        ! Main program
        !=============
        ! Check
        if (size(pe, 1) /= ne) then 
            call gdErrorHandler('ExtractPolygonVertices: ne should equal number of edges')
        end if 
        if (size(pe,2) /= 2) then
            ! Throw error
            call gdErrorHandler('ExtractPolygonVertices: input argument pe should be a ne-by-2 integer array')
    
        end if
        allocate(pv(ne+1))
        pv = 0
        !if (size(pv) /= ne+1) then 
        !    call gdErrorHandler('ExtractPolygonVertices: pv should have dimension ne+1')
        !end if
    
        ! Initialize
        check(:) = .false. 
        k = 1
    
        ! Hedge for the trivial case of one edge
        if (ne == 1) then 
            ! Check 
            if (pe(1,1) == pe(1,2)) then 
                call gdErrorHandler('ExtractPolygonVertices: starting polygon edge vertices are the same - check input')
            else
                pv(1:2) = pe(1,1:2)
            end if 
    
            ! Exit routine
            return 
        end if
    
        ! Get the starting vertex & do check on polygon
        check(1) = pe(1,1) == pe(2,1)
        check(2) = pe(1,1) == pe(2,2)
        check(3) = pe(1,2) == pe(2,1)
        check(4) = pe(1,2) == pe(2,2)
        if ( (check(1) .and. (.not. any(check(2:4))) ) .or. &
             (check(2) .and. (.not. any(check((/1, 3, 4/)))) ) ) then 
            ! Second vertex is starting vertex
            cvind = 2
        else if ( (check(4) .and. (.not. any(check(1:3))) ) .or. &
            (check(3) .and. (.not. any(check((/1, 2, 4/)))) ) ) then 
            ! First vertex is starting vertex
                cvind = 1
        else 
            ! Something wrong with input - throw error
            print *, pe(1,:)
            print *, pe(2,:)
            call gdErrorHandler('ExtractPolygonVertices: something wrong with input - please check input polygon')
        end if
    
        ! Loop
        do while (k < ne)
            ! Get the current vertex, add it to the polygon vertices
            pv(k) = pe(k,cvind)
    
            ! Get the next vertex
            if (cvind == 1) then
                nv = pe(k,2)
            else
                nv = pe(k,1)
            end if 
    
            ! Update cvind
            if (nv == pe(k+1,1)) then
                cvind = 1
            else if (nv == pe(k+1,2)) then
                cvind = 2
            else
                ! Throw error, next edge does not contain the vertex
                call gdErrorHandler('ExtractPolygonVertices: could not find next edge - please check input polygon')
            end if
    
            ! Update counter
            k = k+1
    
        end do
    
        ! Add last vertices
        pv(k) = pe(ne,cvind)
        if (cvind == 1) then
            pv(k+1) = pe(ne,2)
        else
            pv(k+1) = pe(ne,1)
        end if
    
    end subroutine

    ! Vertex and face extractor from unsorted point data
    subroutine ExtractEdgesFromCoordinates(x1, y1, x2, y2, &
        v1, v2, x, y)

        ! Description
        !============
        ! This routine extracts vertex IDs (from 1 to number of vertices
        ! ) for a given set of edge coordinates. Duplicate vertices 
        ! are recognized if the coordinates are the same up to 
        ! their precision (i.e. we use the equals sign to check if 
        ! coordinates are equal). We return v1 and v2 that are the 
        ! vertex indices into x, y coordinates

        ! Declare variables
        !==================
        ! Arguments
        real(R8), dimension(:), intent(in)      :: x1, x2, y1, y2
        integer(I8), dimension(:), allocatable, intent(out) :: v1, v2
        real(R8), allocatable, dimension(:), intent(out)    :: x, y 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: tx, ty, txnew, tynew, &
            thisy, dx
        integer(I8)                             :: si, ei
        integer(I8), allocatable, dimension(:)  :: tv, tvnew, sortind, &
            sortind2, thisind
        
        ! Loop
        integer(I8)                             :: i, j, k, cc  
        
        ! Initialize
        !===========
        ! Sanity checks on dimensions
        if ((size(x1) /= size(y1)) .or. (size(x2) /= size(y2)) .or. &
            (size(x1) /= size(x2)) .or. (size(y1) /= size(y2))) then 
            call gdErrorHandler('ExtractEdgseFromCoordinates: ' // & 
                'incompatible input sizes')
        end if 

        ! Allocate  
        allocate(v1(size(x1)), v2(size(x1)), x(2*size(x1)), y(2*size(x1)))

        ! Initialize
        v1 = 0
        v2 = 0
        x = 0.0_R8
        y = 0.0_R8 

        ! Sort
        !=====
        ! X-coordinates only
        tx = [x1, x2]
        ty = [y1, y2]
        tv = [(k, k = 1, size(tx))]
        allocate(sortind(size(tv)))
        call Sort(tx, ind=sortind, ascend=.true.)
        ty = ty(sortind)
        tv = tv(sortind)

        ! Check which parts contain multiple x-coordinates that are
        ! the same - sort these based on y-value (and adjust sortind)
        dx = tx(2:) - tx(1:size(tx)-1)
        si = 0
        ei = 0
        k = 0
        do while (k < size(dx))
            ! Update counter
            k = k + 1

            ! Check if we need to sort
            if (dx(k) /= 0) then 
                cycle
            end if

            ! Set start
            si = k 

            ! Find end
            ei = findloc(dx(si+1:) /= 0, .true., 1, back=.false.)
            if (ei == 0) then 
                ! End of polygon reached
                ei = size(tx)
            else
                ! Account for starting at si+1
                ei = ei + si
            end if 

            ! Build the index
            thisind = [(j, j = si, ei)]

            ! Sort these coordinates according to y
            allocate(sortind2(size(thisind)))
            thisy = ty(thisind)
            call Sort(thisy, ind=sortind2, ascend=.true.)
            tx(thisind) = tx(thisind(sortind2))
            ty(thisind) = ty(thisind(sortind2))
            sortind(thisind) = sortind(thisind(sortind2))
            deallocate(sortind2)

            ! Update k 
            k = ei - 1
        end do

        ! Check for exactly the same coordinates
        txnew = tx 
        tynew = ty 
        tvnew = tv 
        cc = 1 ! First vertex is always kept
        tvnew(1) = cc 
        do i = 2, size(tx)
            if ((txnew(i) == txnew(i-1)) .and. (tynew(i) == tynew(i-1))) then 
                ! Vertex is the same
                tvnew(i) = tvnew(i-1)
            else
                ! Vertex is different
                cc = cc + 1
                tvnew(i) = cc 
                txnew(cc) = tx(i)
                tynew(cc) = ty(i)
            end if 
        end do 

        ! Construct output
        x = txnew(1:cc)
        y = tynew(1:cc)
        do i = 1, size(tvnew)
            if (sortind(i) > size(x1)) then 
                v2(sortind(i) - size(x1)) = tvnew(i)
            else
                v1(sortind(i)) = tvnew(i)
            end if 
        end do 

    end subroutine

    !------------------------------------------------------------------!
    !                              Numerics                            !
    !------------------------------------------------------------------!

    ! Euclidean distance between two points
    subroutine Distance0D(d, x1, y1, x2, y2) 

        ! Description
        !============
        ! Compute the euclidean distance between two points

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: x1, x2, y1, y2 
        real(R8), intent(out)           :: d

        ! Compute
        !========
        d = sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end subroutine

    ! Euclidean distance between two sets of points (element-wise)
    subroutine Distance1D(d, x1, y1, x2, y2) 

        ! Description
        !============
        ! Compute the euclidean distance between two points

        ! Declare variables
        !==================
        ! Arguments
        real(R8),  intent(in)               ::  x1(:), x2(:), y1(:), y2(:) 
        real(R8), allocatable               :: d(:)

        ! Compute
        !========
        d = sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end subroutine

    ! Continuous index computation
    function ComputeI(x, y, x1, y1, x2, y2) result(frac)

        ! Description
        !============
        ! Compute the relative distance on the edge where a vertex (x, y) lies. It
        ! is assumed that x, y truly lies on the edge, i.e. this routine will give
        ! wrong results if x, y is not on the edge with points (x1, y1), (x2, y2).
        ! Note that also the order of the points matters. 

        ! We hedge for points on one of the two nodes by comparing vertex values
        ! with disttol (this should be conform how the intersections are computed
        ! in SegmentIntersections).

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)  :: x, y, x1, y1, x2, y2
        real(R8), allocatable               :: frac(:)

        ! Auxiliary
        real(R8), allocatable, dimension(:) :: d1, d2, de

        ! Check
        !======
        ! Check for equal distance
        call Distance(d1, x, y, x1, y1)
        call Distance(d2, x, y, x2, y2)
        call Distance(de, x1, y1, x2, y2)

        ! Compute frac
        !=============
        allocate(frac(size(x)))
        where ((d1 < disttol)) frac = 0 ! first point is the same
        where ((d2 < disttol)) frac = 1 ! second point is the same
        where ( .not. (d1 < disttol) .and. .not. (d2 < disttol)) frac = d1/de
        where(frac > 1.0_R8) frac = 1.0_R8
        where(frac < 0.0_R8) frac = 0.0_R8
        !if (any(frac > 1.0_R8) .or. any(frac < 0.0_R8) ) then 
        !    print *, 'roundoff errors detected'
        !end if 

    end function 

    ! Intersection between two lines with points x11, x12, x21, x22 
    subroutine LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! Compute the intersection between two lines. Hereto, we solve
        ! the following system for x, y:
        !
        !       (y12 - y11)(x - x11) = (x12 - x11)(y - y11)
        !       (y22 - y21)(x - x21) = (x22 - x21)(y - y21)
        !
        ! The solution can be computed analytically in this case 
        ! (provided some conditions on the differences, which is checked
        ! for with disttol) 

        ! If the intersection does not exist (e.g. parallel lines up to
        ! precision disttol),  x and y are NaN valued. If the lines are
        ! collinear, we return inf. The latter is checked by computing 
        ! the normal distance between the two lines (which, again, is 
        ! checked by computing the normal to one of the lines and 
        ! projecting the vector between two points of both lines onto 
        ! this normal and checking the length). If this distance is 
        ! lower than the distance tolerance, the lines are said to be 
        ! collinear. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: det, r1, r2, nan, dx1, dx2, dy1, &
            dy2, dist, nx, ny, nn, vx, vy, inf 

        ! Loop

        ! Initialize
        !===========
        ! Compute
        dx1 = x12 - x11 
        dy1 = y12 - y11 
        dx2 = x22 - x21 
        dy2 = y22 - y21 

        ! Some checks
        if ( (abs(dx1) < disttol) .and. (abs(dy1) < disttol) ) then 
            ! Print warning
            call PolygonWarningHandler('LineIntersections: ' &
                // 'first line is a point up to distance precision, ' &
                // 'results may be inaccurate. Consider rescaling the' &
                // ' coordinates')
        end if 
        if ( (abs(dx2) < disttol) .and. (abs(dy2) < disttol) ) then 
            ! Print warning
            call PolygonWarningHandler('LineIntersections: ' &
                // 'second line is a point up to distance precision, ' &
                // 'results may be inaccurate. Consider rescaling the' &
                // ' coordinates')
        end if 

        ! Check determinant/slopes
        r1 = dx1/dy1
        r2 = dx2/dy2 

        det = r1 - r2 
        if ( (abs(det) < macheps) .or. (abs(dy1) < macheps .and. abs(dy2) < macheps)) then 
            ! Parallel or collinear lines - need to check collinearity
            
            ! Compute normal
            nx = -(y11 - y12)
            ny = (x11 - x12)
            nn = sqrt(nx**2 + ny**2)

            ! Compute vector between lines
            vx = (x11 - x21)
            vy = (y11 - y21)

            ! Compute the distance
            dist = abs( vx*nx/nn + vy*ny/nn )

            ! Check 
            if (dist < disttol) then 
                ! collinear lines, return inf
                x = IEEE_VALUE(inf, IEEE_positive_inf)
                y = x 
            else 
                ! Parallel lines, return nan
                x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                y = x 
            end if 
            return 
        end if 

        ! Compute intersection
        !=====================
        ! Hedge for small dx when computing slope 
        if (abs(dx1) > abs(dx2)) then 
            if (abs(dx1) > disttol) then 
                r1 = dy1/dx1 
                if (abs(dx2) > disttol) then 
                    ! Two non-parallel, non-vertical and non-horizontal lines
                    r2 = dy2/dx2 
                    x = (r1*x11 - r2*x21 -y11 + y21)/(r1 - r2)
                    y = r1*(x - x11) + y11
                else 
                    ! Second line is vertical line, first one is non-vertical
                    x = x21
                    y = r1*(x - x11) + y11
                end if 

            else
                ! Both lines are parallel - should've been captured before actually
                ! Compute normal
                nx = -(y11 - y12)
                ny = (x11 - x12)
                nn = sqrt(nx**2 + ny**2)

                ! Compute vector between lines
                vx = (x11 - x21)
                vy = (y11 - y21)

                ! Compute the distance
                dist = abs( vx*nx/nn + vy*ny/nn )

                ! Check 
                if (dist < disttol) then 
                    ! collinear lines, return inf
                    x = IEEE_VALUE(inf, IEEE_positive_inf)
                    y = x 
                else 
                    ! Parallel lines, return nan
                    x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                    y = x 
                end if 
                return 
            end if 
        else
            if (abs(dx2) > disttol) then 
                r2 = dy2/dx2
                if (abs(dx1) > disttol) then 
                    ! Two non-parallel, non-vertical and non-horizontal lines
                    r1 = dy1/dx1 
                    x = (r1*x11 - r2*x21 -y11 + y21)/(r1 - r2)
                    y = r2*(x - x21) + y21
                else 
                    ! First line is vertical line, second one is non-vertical
                    x = x11
                    y = r2*(x - x21) + y21
                end if 

            else
                ! Both lines are parallel - should've been captured before actually
                ! Compute normal
                nx = -(y11 - y12)
                ny = (x11 - x12)
                nn = sqrt(nx**2 + ny**2)

                ! Compute vector between lines
                vx = (x11 - x21)
                vy = (y11 - y21)

                ! Compute the distance
                dist = abs( vx*nx/nn + vy*ny/nn )

                ! Check 
                if (dist < disttol) then 
                    ! collinear lines, return inf
                    x = IEEE_VALUE(inf, IEEE_positive_inf)
                    y = x 
                else 
                    ! Parallel lines, return nan
                    x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                    y = x 
                end if 
                return 
            end if 
        end if 


    end subroutine

    ! Intersection between segment and line 
    subroutine SegmentLineIntersecions(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! This routine computes the intersection between a line segment
        ! given by the points (x11, y11) and (x12, y12) and the line
        ! given by the points (x21, y21), (x22, y22). First, it is 
        ! checked whether there is an intersection between the two lines
        ! formed by this point. If there is an intersection, then it is
        ! checked whether this intersection lies on the line segment. 
        ! The latter is done by checking the scalar product between the
        ! vectors formed by the intersection point and the segment's 
        ! points (this should be negative if the point lies on the 
        ! segment). We hedge for machine precision effects by checking
        ! if the point lies within macheps to one of the nodes. 
        
        ! If no intersection is found, NaNs are returned for x and y

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: nan, d1, d2, dotprod

        ! Compute intersection
        !=====================
        call LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
            x22, y22)

        ! Return if no intersection is found - x and y will be NaN 
        ! already
        if (isnan(x)) then 
            return 
        end if 

        ! Check if intersection is on segment
        dotprod = (x - x11)*(x - x12) + (y - y11)*(y - y12)
        call Distance(d1, x, y, x11, y11) 
        call Distance(d2, x, y, x12, y12)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            return 
        end if 
        if (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
        end if 

    end subroutine

    ! Intersection between two segments
    subroutine SegmentIntersections(x, y, x11, y11, x12, y12, x21, y21, &
        x22, y22)

        ! Description
        !============
        ! This routine computes the intersection between two line segments
        ! given by the points (x11, y11) and (x12, y12) and 
        ! (x21, y21), (x22, y22). First, it is 
        ! checked whether there is an intersection between the two lines
        ! formed by these points. If there is an intersection, then it is
        ! checked whether this intersection lies on both line segments. 
        ! The latter is done by checking the scalar product between the
        ! vectors formed by the intersection point and the segment's 
        ! points (this should be negative if the point lies on the 
        ! segment). We hedge for machine precision effects by checking
        ! if the point lies within macheps to one of the nodes. 
        
        ! If no intersection is found, NaNs are returned for x and y. 
        ! For segments that are collinear and overlap, there are in 
        ! principle an infinite amount of intersections. However, for 
        ! these we only return a single intersection located in the 
        ! middle of the overlapping part of the segment. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: nan, d1, d2, dotprod, d11, d12, &
            d21, d22, dp11, dp12, dp21, dp22, vx1121, vx1122, vx1221, vx1222, &
            vy1121, vy1122, vy1221, vy1222, x2, y2
        logical                 :: v11on2, v12on2, v21on1, v22on1, hasoverlap

        ! Compute intersection
        !=====================
        ! Check if edges overlap, if so: compute intersections
        hasoverlap = CheckEdgeOverlap(x11, y11, x12, y12, x21, y21, x22, y22)
        if (hasoverlap) then 
            call LineIntersections(x, y, x11, y11, x12, y12, x21, y21, &
                x22, y22)
        else 
            ! No overlap -> no intersections
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x
        end if 

        ! First, check if nodes coincide - overwrite in that case (to hedge for
        ! ill-conditioning effects in line intersection routine)
        call CheckSegmentCoincidence(x2, y2, x11, y11, x12, y12, x21, y21, &
            x22, y22)
        if (.not. isnan(x2)) then 
            x = x2
            y = y2
        end if

        ! Return if no intersection is found - x and y will be NaN 
        ! already
        if (isnan(x)) then 
            return 
        end if 

        ! If the result is inf, the segments are collinear. Then we need
        ! to check if the edges overlap. If so, the intersection is 
        ! located at the middle of the overlapping section (if the 
        ! intersection is in a point, then this point should be returned)
        if (.not. ieee_is_finite(x)) then 
            ! Distances between segment points
            call Distance(d11, x11, y11, x21, y21) 
            call Distance(d12, x11, y11, x22, y22)
            call Distance(d21, x12, y12, x21, y21) 
            call Distance(d22, x12, y12, x22, y22)

            ! Check 
            if ( ((d11 < disttol) .and. (d22 < disttol)) &
                 .or. ((d12 < disttol) .and. (d21 < disttol)) ) then 
                ! Edges are the same up to distance tolerance, take
                ! average - should actually already be captured before
                x = 0.5*(x11 + x12) 
                y = 0.5*(y11 + y12) 
            else
                ! Compute vectors between points 
                vx1121 = x21 - x11 
                vx1122 = x22 - x11 
                vx1221 = x21 - x12 
                vx1222 = x22 - x12 

                vy1121 = y21 - y11 
                vy1122 = y22 - y11 
                vy1221 = y21 - y12 
                vy1222 = y22 - y12 

                ! Compute scalar products
                dp11 = vx1121*vx1122 + vy1121*vy1122 
                dp12 = vx1221*vx1222 + vy1221*vy1222 

                dp21 = vx1221*vx1121 + vy1221*vy1121 
                dp22 = vx1222*vx1122 + vy1222*vy1122

                ! Checks
                if (d11 < disttol) then 
                    ! First points coincide 
                    if (dp12 < 0) then 
                        ! Second point of first segment on second 
                        ! segment, take average of first segment
                        x = 0.5*(x11 + x12) 
                        y = 0.5*(y11 + y12) 
                    elseif (dp22 < 0) then 
                        ! Second point of second segment on first 
                        ! segment, take average of second segment
                        x = 0.5*(x21 + x22) 
                        y = 0.5*(y21 + y22)
                    else    
                        ! Intersection in single point 
                        x = x11 
                        y = y11 
                    end if 
                elseif (d12 < disttol) then 
                    ! First point of first segment and second point of 
                    ! second segment coincide
                    if (dp12 < 0) then 
                        ! Second point of first segment on second 
                        ! segment, take average of first segment
                        x = 0.5*(x11 + x12) 
                        y = 0.5*(y11 + y12) 
                    elseif (dp21 < 0) then 
                        ! First point of second segment on first 
                        ! segment, take average of second segment
                        x = 0.5*(x21 + x22) 
                        y = 0.5*(y21 + y22)
                    else    
                        ! Intersection in single point
                        x = x11 
                        y = y11
                    end if 
                elseif (d21 < disttol) then 
                    ! Second point of first segment and first point of 
                    ! second segment coincide
                    if (dp11 < 0) then 
                        ! First point of first segment on second 
                        ! segment, take average of first segment
                        x = 0.5*(x11 + x12) 
                        y = 0.5*(y11 + y12) 
                    elseif (dp22 < 0) then 
                        ! Second point of second segment on first 
                        ! segment, take average of second segment
                        x = 0.5*(x21 + x22) 
                        y = 0.5*(y21 + y22)
                    else    
                        ! Intersection in single point
                        x = x12 
                        y = y12
                    end if 
                elseif (d22 < disttol) then 
                    ! Second point of first segment and first point of 
                    ! second segment coincide
                    if (dp11 < 0) then 
                        ! First point of first segment on second 
                        ! segment, take average of first segment
                        x = 0.5*(x11 + x12) 
                        y = 0.5*(y11 + y12) 
                    elseif (dp21 < 0) then 
                        ! First point of second segment on first 
                        ! segment, take average of second segment
                        x = 0.5*(x21 + x22) 
                        y = 0.5*(y21 + y22)
                    else    
                        ! Intersection in single point
                        x = x12 
                        y = y12
                    end if 
                else 
                    ! No coincidence of nodes, can simply check vector 
                    ! products
                    if (vx1121*vx1122 + vy1121*vy1122 < 0) then 
                        ! Node 11 lies on segment 2
                        v11on2 = .true.
                    end if 
                    if (vx1221*vx1222 + vy1221*vy1222 < 0) then 
                        ! Node 12 lies on segment 2
                        v12on2 = .true.
                    end if 
                    if (vx1121*vx1221 + vy1121*vy1221 < 0) then 
                        ! Node 21 lies on segment 1
                        v21on1 = .true.
                    end if 
                    if (vx1122*vx1221 + vy1122*vy1221 < 0) then 
                        ! Node 22 lies on segment 1
                        v22on1 = .true.
                    end if 

                    ! Check cases
                    if (v11on2 .and. v12on2) then 
                        ! Segment 1 lies within segment 2
                        x = 0.5*(x11 + x12)
                        y = 0.5*(y11 + y12)
                    elseif (v21on1 .and. v22on1) then 
                        ! Segment 2 lies within segment 1
                        x = 0.5*(x21 + x22)
                        y = 0.5*(y21 + y22)
                    elseif (v11on2 .and. v21on1) then 
                        x = 0.5*(x11 + x21)
                        y = 0.5*(y11 + y21)
                    elseif (v11on2 .and. v22on1) then 
                        x = 0.5*(x11 + x22)
                        y = 0.5*(y11 + y22)
                    elseif (v12on2 .and. v21on1) then 
                        x = 0.5*(x12 + x21)
                        y = 0.5*(y12 + y21)
                    elseif (v12on2 .and. v22on1) then 
                        x = 0.5*(x12 + x22)
                        y = 0.5*(y12 + y22)
                    else 
                        ! Edges don't seem to be overlapping, shouldn't
                        ! happen at this point. Throw error
                        call PolygonErrorHandler('Edges that should overlap do not seem to overlap')
                    end if 
                end if 
            end if 
            return
        end if

        ! Check if intersection is on segment 1
        dotprod = (x - x11)*(x - x12) + (y - y11)*(y - y12)
        call Distance(d1, x, y, x11, y11) 
        call Distance(d2, x, y, x12, y12)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            ! Lies on vertex here, still need to check other segment
        elseif (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
            return
        end if 
    
        ! Check if intersection is on segment 2
        dotprod = (x - x21)*(x - x22) + (y - y21)*(y - y22)
        call Distance(d1, x, y, x21, y21) 
        call Distance(d2, x, y, x22, y22)
        if ( (d1 < disttol) .or. (d2 < disttol) ) then 
            ! Lies on vertex here
        elseif (dotprod > 0) then 
            x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            y = x 
            return 
        end if 

    end subroutine

    ! Intersections between a segment and a polygon
    subroutine SegmentPolygonIntersections(polygon, x1, y1, x2, y2, &
        x, y, s)

        ! Description
        !============
        ! Compute the intersections between a segment and a polygon. The
        ! segment must be given by two points, (x1, y1), (x2, y2). 
        ! The output in x, y are the intersection coordinates. s 
        ! contains a list (integer) of polygon segments where 
        ! intersections where found. 
        
        ! Algorithm
        !==========
        ! The main idea is to loop over all segments of the polygon and
        ! check for polygon edge whether it has an intersection with 
        ! the segment. To compute intersections, the 
        ! SegmentIntersections routine is used. To reduce computational 
        ! cost, a simple check is made whether the encompassing boxes of
        ! two edges overlap or not. If they don't, there can be no 
        ! intersection. 
        
        ! Notes
        !======
        ! Note 1: we now hedge for duplicate intersections that happen 
        ! exactly in one of the nodes of the polygon, and which appears
        ! twice due to the fact that the node belongs to two edges, for
        ! which intersections are sought. Note that actual multiple 
        ! intersections (same coordinates, but different segment 
        ! indices) are not removed, since these are valid intersections. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT), intent(in)           :: polygon 
        real(R8), intent(in)                    :: x1, y1, x2, y2 
        real(R8), allocatable, intent(out)      :: x(:), y(:) 
        integer(I8), allocatable, intent(out)   :: s(:)

        ! Auxiliary
        integer(I8)                         :: counter  
        real(R8)                            :: xi, yi, xe1, ye1, xe2, ye2, &
            d

        integer(I8), allocatable            :: temps(:), sortindex(:)
        real(R8), allocatable               :: tempx(:), tempy(:), &
            tempsr(:)
        logical, allocatable                :: keepind(:)

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then 
            deallocate(x) 
        end if    
        if (allocated(y)) then 
            deallocate(y) 
        end if        

        ! Unpack polygon
        associate( &
            ne          => polygon%ne,      & 
            edges       => polygon%edges,   & 
            xp          => polygon%x,       &
            yp          => polygon%y)

        ! Initialize intersection counter
        counter = 0

        ! Allocate temporary arrays
        allocate(tempx(ne), tempy(ne), temps(ne), tempsr(ne)) ! maximum ne intersections, to be trimmed later

        ! Loop   
        do i = 1, ne 
            ! Get coordinates of next polygon edge
            xe1 = xp(edges(i, 1))
            ye1 = yp(edges(i, 1))
            xe2 = xp(edges(i, 2))
            ye2 = yp(edges(i, 2))

            ! Check boxes
            if (CheckEdgeOverlap(xe1, ye1, xe2, ye2, x1, y1, x2, y2)) then 
                ! Edges overlap, compute intersection
                call SegmentIntersections(xi, yi, xe1, ye1, xe2, ye2, &
                    x1, y1, x2, y2)
                
                ! If intersection is found, add it
                if (.not. isnan(xi)) then 
                    ! Update counter
                    counter = counter + 1
                    tempx(counter) = xi 
                    tempy(counter) = yi 
                    temps(counter) = i
                end if 
            end if
        end do  

        ! Allocate and attribute
        allocate(x(counter), y(counter), s(counter))
        x = tempx(1:counter) 
        y = tempy(1:counter)  
        s = temps(1:counter)

        ! Hedge for duplicates
        !=====================
        ! Since intersections should be sorted by default, we can 
        ! simply loop and check

        ! Sort intersections along s - may not be sorted anymore due to 
        ! parallellization
        allocate(sortindex(counter))
        call Sort(s, ind=sortindex)
        x = x(sortindex)
        y = y(sortindex)

        ! Hedge for duplicates
        allocate(keepind(counter))
        keepind = .true.
        do i = 1, counter-1
            ! Check for segment index
            if ((s(i+1) - s(i)) == 1) then 
                ! Check for same intersection
                call Distance(d, x(i), y(i), x(i+1), y(i+1))
                if (d <= disttol) then 
                    ! Delete
                    keepind(i+1) = .false. 
                end if 
            end if 
        end do 

        ! Delete
        x = pack(x, keepind)
        y = pack(y, keepind)
        s = pack(s, keepind)

        ! Housekeeping
        !=============
        end associate

        deallocate(tempx, tempy)



    end subroutine

    subroutine SegmentSimplePolygonIntersections(xp, yp, x1, y1, x2, y2, &
        x, y, s)

        ! Description
        !============
        ! Compute the intersections between a segment and a polygon. The
        ! segment must be given by two points, (x1, y1), (x2, y2). 
        ! The output in x, y are the intersection coordinates. s 
        ! contains a list (integer) of polygon segments where 
        ! intersections where found. 
        
        ! Algorithm
        !==========
        ! The main idea is to loop over all segments of the polygon and
        ! check for polygon edge whether it has an intersection with 
        ! the segment. To compute intersections, the 
        ! SegmentIntersections routine is used. To reduce computational 
        ! cost, a simple check is made whether the encompassing boxes of
        ! two edges overlap or not. If they don't, there can be no 
        ! intersection. 
        
        ! Notes
        !======
        ! Note 1: we now hedge for duplicate intersections that happen 
        ! exactly in one of the nodes of the polygon, and which appears
        ! twice due to the fact that the node belongs to two edges, for
        ! which intersections are sought. Note that actual multiple 
        ! intersections (same coordinates, but different segment 
        ! indices) are not removed, since these are valid intersections. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)      :: xp, yp
        real(R8), intent(in)                    :: x1, y1, x2, y2 
        real(R8), allocatable, intent(out)      :: x(:), y(:) 
        integer(I8), allocatable, intent(out)   :: s(:)

        ! Auxiliary
        integer(I8)                         :: counter, ne
        real(R8)                            :: xi, yi, xe1, ye1, xe2, ye2, &
            d

        integer(I8), allocatable            :: temps(:), sortindex(:)
        real(R8), allocatable               :: tempx(:), tempy(:), &
            tempsr(:)
        logical, allocatable                :: keepind(:)

        ! Loop
        integer(I8)                         :: i 

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then 
            deallocate(x) 
        end if    
        if (allocated(y)) then 
            deallocate(y) 
        end if        

        ! Initialize intersection counter
        counter = 0

        ! Precompute overlapping edges
        ne = size(xp)-1

        ! Allocate temporary arrays
        allocate(tempx(ne), tempy(ne), temps(ne), tempsr(ne)) ! maximum ne intersections, to be trimmed later

        ! Loop   
        do i = 1, ne
            ! Get coordinates of next polygon edge
            xe1 = xp(i)
            ye1 = yp(i)
            xe2 = xp(i+1)
            ye2 = yp(i+1)

            ! Check boxes
            if (CheckEdgeOverlap(xe1, ye1, xe2, ye2, x1, y1, x2, y2)) then 
                ! Edges overlap, compute intersection
                call SegmentIntersections(xi, yi, xe1, ye1, xe2, ye2, &
                    x1, y1, x2, y2)
                
                ! If intersection is found, add it
                if (.not. isnan(xi)) then 
                    ! Update counter
                    counter = counter + 1
                    tempx(counter) = xi 
                    tempy(counter) = yi 
                    temps(counter) = i
                end if 
            end if
        end do  

        ! Allocate and attribute
        allocate(x(counter), y(counter), s(counter))
        x = tempx(1:counter)
        y = tempy(1:counter)
        s = temps(1:counter)

        ! Hedge for duplicates
        !=====================
        ! Since intersections should be sorted by default, we can 
        ! simply loop and check

        ! Sort intersections along s - may not be sorted anymore due to 
        ! parallellization
        allocate(sortindex(counter))
        call Sort(s, ind=sortindex)
        x = x(sortindex)
        y = y(sortindex)

        ! Hedge for duplicates
        allocate(keepind(counter))
        keepind = .true.
        do i = 1, counter-1
            ! Check for segment index
            if ((s(i+1) - s(i)) == 1) then 
                ! Check for same intersection
                call Distance(d, x(i), y(i), x(i+1), y(i+1))
                if (d <= disttol) then 
                    ! Delete
                    keepind(i+1) = .false. 
                end if 
            end if 
        end do 

        ! Delete
        x = pack(x, keepind)
        y = pack(y, keepind)
        s = pack(s, keepind)

        ! Housekeeping
        !=============
        deallocate(tempx, tempy)



    end subroutine

    ! Intersections of polygon with itself
    subroutine PolygonSelfIntersections(polygon, x, y, s1, s2)

        ! Description
        !============
        ! This routine computes the intersections of a polygon with 
        ! itself. Hereto, we loop over all polygon edges and compute
        ! intersections with all next edges. Note that intersections 
        ! between neighbouring edges are not considered to be self-
        ! intersections! Also for closed polygons, the 'intersection'
        ! of the last and first edge is not accounted for. It is assumed
        ! that all polygon fields are up to date. 

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT), intent(in)           :: polygon 
        real(R8), allocatable, intent(out)      :: x(:), y(:)
        integer(I8), allocatable, intent(out)   :: s1(:), s2(:)

        ! Auxiliary
        integer(I8)                             :: counter
        real(R8)                                :: xei1, yei1, xei2, &
            yei2, xej1, xej2, yej1, yej2, xi, yi
        real(R8), allocatable                   :: tempx(:), tempy(:)
        integer(I8), allocatable                :: temps1(:), &
            temps2(:)

        ! Loop
        integer(I8)                             :: i, j 

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if
        if (allocated(s1)) then 
            deallocate(s1) 
        end if
        if (allocated(s2)) then 
            deallocate(s2) 
        end if

        ! Associate
        associate(&
            xp          => polygon%x,        &
            yp          => polygon%y,        &
            edges       => polygon%edges,    &
            ne          => polygon%ne)

        ! Initialize
        counter     = 0 ! intersection counter 

        ! Check for trivial inputs
        if (ne <= 2) then 
            ! No intersections possible
            allocate(x(0), y(0), s1(0), s2(0))
            return 
        end if

        ! Allocate (too large, trim later)
        allocate(tempx(ne), tempy(ne), temps1(ne), temps2(ne))

        ! Compute intersections
        !======================
        ! Loop over the edges of the polygon
        do i = 2, ne-2 
            ! Get coordinates of this edge
            xei1 = xp(edges(i, 1))
            yei1 = yp(edges(i, 1))
            xei2 = xp(edges(i, 2))
            yei2 = yp(edges(i, 2))

            ! Loop over remaining polygon segments, but skip 
            ! neighbouring edges 
            do j = i+3, ne 
                ! Get coordinates of this edge
                xej1 = xp(edges(j, 1))
                yej1 = yp(edges(j, 1))
                xej2 = xp(edges(j, 2))
                yej2 = yp(edges(j, 2))
                
                ! Compute intersections
                call SegmentIntersections(xi, yi, xei1, yei1, xei2, yei2, &
                    xej1, yej1, xej2, yej2) 
                
                ! Check if there is an intersection 
                if (.not. isnan(xi)) then  
                    ! Intersection found, add
                    counter = counter + 1
                    if (counter > size(tempx)) then 
                        ! Extend
                        tempx = [tempx, 0*tempx]
                        tempy = [tempy, 0*tempy]
                        temps1 = [temps1, 0*temps1]
                        temps2 = [temps2, 0*temps2]
                    end if 
                    tempx(counter) = xi 
                    tempy(counter) = yi 
                    temps1(counter) = i 
                    temps2(counter) = j 
                end if 
            end do 
        end do 

        ! If the polygon is not closed, check the last edge and first 
        ! edge
        if (.not. polygon%isclosed) then 
            i = 1
            j = ne 
            xei1 = xp(edges(i, 1))
            yei1 = yp(edges(i, 1))
            xei2 = xp(edges(i, 2))
            yei2 = yp(edges(i, 2))
            xej1 = xp(edges(j, 1))
            yej1 = yp(edges(j, 1))
            xej2 = xp(edges(j, 2))
            yej2 = yp(edges(j, 2))
            call SegmentIntersections(xi, yi, xei1, yei1, xei2, yei2, &
                xej1, yej1, xej2, yej2) 
            
            
            ! Check if there is an intersection 
            if (.not. isnan(xi)) then  
                ! Intersection found, add
                counter = counter + 1
                if (counter > size(tempx)) then 
                    ! Extend
                    tempx = [tempx, 0*tempx]
                    tempy = [tempy, 0*tempy]
                    temps1 = [temps1, 0*temps1]
                    temps2 = [temps2, 0*temps2]
                end if 
                tempx(counter) = xi 
                tempy(counter) = yi 
                temps1(counter) = i 
                temps2(counter) = j 
            end if

            if (ne > 3) then 
                j = 3 
                xei1 = xp(edges(i, 1))
                yei1 = yp(edges(i, 1))
                xei2 = xp(edges(i, 2))
                yei2 = yp(edges(i, 2))
                xej1 = xp(edges(j, 1))
                yej1 = yp(edges(j, 1))
                xej2 = xp(edges(j, 2))
                yej2 = yp(edges(j, 2))
                call SegmentIntersections(xi, yi, xei1, yei1, xei2, yei2, &
                    xej1, yej1, xej2, yej2) 
                
                
                ! Check if there is an intersection 
                if (.not. isnan(xi)) then  
                    ! Intersection found, add
                    counter = counter + 1
                    if (counter > size(tempx)) then 
                        ! Extend
                        tempx = [tempx, 0*tempx]
                        tempy = [tempy, 0*tempy]
                        temps1 = [temps1, 0*temps1]
                        temps2 = [temps2, 0*temps2]
                    end if 
                    tempx(counter) = xi 
                    tempy(counter) = yi 
                    temps1(counter) = i 
                    temps2(counter) = j 
                end if  
            end if 
        end if 
            

        ! Add to output
        allocate(x(counter), y(counter), s1(counter), s2(counter)) 
        x = tempx(1:counter) 
        y = tempy(1:counter) 
        s1 = temps1(1:counter) 
        s2 = temps2(1:counter) 

        ! Housekeeping
        !=============
        end associate 

        deallocate(tempx, tempy, temps1, temps2)

    end subroutine

    ! Intersections between two polygons
    subroutine PolygonIntersections(p1, p2, x, y, s1, s2, s1r, s2r)

        ! Description
        !============
        ! This routine computes the intersections between two polygons.
        ! It returns the coordinates of these intersections in the 
        ! x, y arrays and the edge numbers in s1 and s2 for the first 
        ! and second polygon, resp. 

        ! Important: don't use this routine to compute
        ! self-intersections! Do this with the dedicated 
        ! PolygonSelfIntersections routine... 

        ! Algorithm
        !==========
        ! We simply loop over all edges of p2 and call 
        ! SegmentPolygonIntersections of p1 to compute the 
        ! intersections. Since the number of intersections is a priori
        ! unknown, and since the maximal amount of intersections may 
        ! be very large (but usually very small), we dynamically grow 
        ! the intersection storage arrays while computing intersections.

        ! Notes
        !======
        ! Note 1: we now hedge for duplicate intersections that happen 
        ! exactly in one of the nodes of the polygon, and which appears
        ! twice due to the fact that the node belongs to two edges, for
        ! which intersections are sought. Note that actual multiple 
        ! intersections (same coordinates, but different segment 
        ! indices) are not removed, since these are valid intersections. 
 
        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT), intent(in)           :: p1, p2
        real(R8), allocatable, intent(out)      :: x(:), y(:)
        integer(I8), allocatable, intent(out)   :: s1(:), s2(:)
        real(R8), allocatable, intent(out), optional    :: s1r(:), s2r(:)

        ! Auxiliary
        integer(I8)                             :: ni, counter, sz, &
            szold, szmult 
        real(R8)                                :: xe1, ye1, xe2, ye2, &
            d
        real(R8), allocatable                   :: tempx(:), tempy(:), &
            xi(:), yi(:), temps1r(:), temps2r(:)
        integer(I8), allocatable                :: temps1(:), &
            temps2(:), si(:), sortindex(:)
        logical, allocatable                    :: keepind(:)

        ! Loop
        integer(I8)                             :: i 

        ! Memory mgmt
        integer(I8), allocatable                :: mgmti(:, :)
        real(R8), allocatable                   :: mgmtr(:, :)

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if
        if (allocated(s1)) then 
            deallocate(s1) 
        end if
        if (allocated(s2)) then 
            deallocate(s2) 
        end if
        if (present(s1r)) then 
            if (allocated(s1r)) then 
                deallocate(s1r)
            end if
        end if 
        if (present(s2r)) then 
            if (allocated(s2r)) then 
                deallocate(s2r)
            end if 
        end if 

        ! Check if either both s1r, s2r are present or not present
        if (((.not. present(s1r)) .and. present(s1r)) &
            .or. (present(s1r) .and. .not. present(s2r))) then 
            call gdErrorHandler('PolygonIntersections: s1r and s2r should either be ' // &
                'both present or not present, one of the two is not supported')
        end if 

        ! Associate
        associate(&
            xp          => p2%x,        &
            yp          => p2%y,        &
            edges       => p2%edges,    &
            ne          => p2%ne)

        ! Initialize
        counter     = 0 ! intersection counter 
        ni          = 0
        szold       = 0
        sz          = 2 ! initial size of intersection array
        szmult      = 2 ! size multiplier 

        ! Allocate
        allocate(tempx(sz), tempy(sz), temps1(sz), temps2(sz), &
            temps1r(sz), temps2r(sz))

        ! Compute intersections
        !======================
        ! Loop over p2
        do i = 1, ne 
            ! Get coordinates of next polygon edge
            xe1 = xp(edges(i, 1))
            ye1 = yp(edges(i, 1))
            xe2 = xp(edges(i, 2))
            ye2 = yp(edges(i, 2))

            ! Compute intersections (duplicates of p1 with segment are 
            ! already removed in this routine)
            call SegmentPolygonIntersections(p1, xe1, ye1, xe2, ye2, &
                xi, yi, si)

            ! Check if intersections were found
            ni = size(xi)
            if (ni > 0) then 
                ! Memory MGMT
                if (counter + ni > sz) then 
                    ! Store old size
                    szold = sz

                    ! Store old values
                    allocate(mgmti(szold, 2), mgmtr(szold, 2))
                    mgmti(:, 1) = temps1 
                    mgmti(:, 2) = temps2   
                    mgmtr(:, 1) = tempx 
                    mgmtr(:, 2) = tempy

                    ! Adjust size
                    do while (sz < counter+ni)
                        sz = sz*szmult 
                    end do

                    ! Reallocate
                    deallocate(tempx, tempy, temps1, temps2, temps1r, temps2r) 
                    allocate(tempx(sz), tempy(sz), temps1(sz), temps2(sz), &
                        temps1r(sz), temps2r(sz))

                    ! Add
                    tempx(1:szold) = mgmtr(:, 1) 
                    tempy(1:szold) = mgmtr(:, 2)
                    temps1(1:szold) = mgmti(:, 1)
                    temps2(1:szold) = mgmti(:, 2)

                    ! Deallocate mgmt arrays
                    deallocate(mgmti, mgmtr)
                end if 

                ! Add intersections
                tempx(counter+1:counter+ni) = xi 
                tempy(counter+1:counter+ni) = yi
                temps1(counter+1:counter+ni) = si 
                temps2(counter+1:counter+ni) = i

                ! Update counter
                counter = counter + ni
            end if 
        end do  

        ! Add to output
        allocate(x(counter), y(counter), s1(counter), s2(counter)) 
        x = tempx(1:counter) 
        y = tempy(1:counter) 
        s1 = temps1(1:counter) 
        s2 = temps2(1:counter) 

        ! Hedge for duplicates
        !=====================
        ! We only need to check s2 since duplicates of s1 should have
        ! already been removed before

        ! Sort intersections along s2 - may not be sorted anymore due to 
        ! parallellization
        allocate(sortindex(counter))
        call Sort(s2, ind=sortindex)
        x = x(sortindex)
        y = y(sortindex)
        s1 = s1(sortindex)

        ! Hedge for duplicates
        allocate(keepind(counter))
        keepind = .true.
        do i = 1, counter-1
            ! Check for segment index
            if ((s2(i+1) - s2(i)) == 1) then 
                ! Check for same intersection
                call Distance(d, x(i), y(i), x(i+1), y(i+1))
                if (d <= disttol) then 
                    ! Delete
                    keepind(i+1) = .false. 
                end if 
            end if 
        end do 

        ! Delete
        x = pack(x, keepind)
        y = pack(y, keepind)
        s1 = pack(s1, keepind)
        s2 = pack(s2, keepind)
        counter = count(keepind)

        ! Compute true intersection locations
        !====================================
        if (present(s1r) .and. present(s2r)) then 
            ! Compute the continuous intersection index (0: first point
            ! of polygon, ne+1: last point of polygon) - note: we need
            ! to use the vert array here instead of the edges array, since
            ! the latter is not necessarily sorted!
            allocate(s1r(counter), s2r(counter))

            ! First polygon index
            !s1r = ComputeI(x, y, p1%x(p1%edges(s1, 1)), p1%y(p1%edges(s1, 1)), &
            !    p1%x(p1%edges(s1, 2)), p1%y(p1%edges(s1, 2))) + s1 - 1 
            s1r = ComputeI(x, y, p1%x(p1%vert(s1)), p1%y(p1%vert(s1)), &
                p1%x(p1%vert(s1+1)), p1%y(p1%vert(s1+1))) + s1 - 1 

            ! Second polygon index
            !s2r = ComputeI(x, y, p2%x(p2%edges(s2, 1)), p2%y(p2%edges(s2, 1)), &
            !    p2%x(p2%edges(s2, 2)), p2%y(p2%edges(s2, 2))) + s2 - 1
            s2r = ComputeI(x, y, p2%x(p2%vert(s2)), p2%y(p2%vert(s2)), &
                p2%x(p2%vert(s2+1)), p2%y(p2%vert(s2+1))) + s2 - 1

        end if 

        ! Housekeeping
        !=============
        end associate 

        deallocate(tempx, tempy, temps1, temps2)


    end subroutine

    ! Intersections between two simple polygons (given by coordinates only)
    subroutine SimplePolygonIntersections(x1, y1, x2, y2, x, y, s1, s2, s1r, s2r)

        ! Description
        !============
        ! This routine computes the intersections between two polygons.
        ! It returns the coordinates of these intersections in the 
        ! x, y arrays and the edge numbers in s1 and s2 for the first 
        ! and second polygon, resp. 

        ! Important: don't use this routine to compute
        ! self-intersections! Do this with the dedicated 
        ! PolygonSelfIntersections routine... 

        ! Algorithm
        !==========
        ! We simply loop over all edges of p2 and call 
        ! SegmentPolygonIntersections of p1 to compute the 
        ! intersections. Since the number of intersections is a priori
        ! unknown, and since the maximal amount of intersections may 
        ! be very large (but usually very small), we dynamically grow 
        ! the intersection storage arrays while computing intersections.

        ! Notes
        !======
        ! Note 1: we now hedge for duplicate intersections that happen 
        ! exactly in one of the nodes of the polygon, and which appears
        ! twice due to the fact that the node belongs to two edges, for
        ! which intersections are sought. Note that actual multiple 
        ! intersections (same coordinates, but different segment 
        ! indices) are not removed, since these are valid intersections. 
 
        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)      :: x1, y1, x2, y2
        real(R8), allocatable, intent(out)      :: x(:), y(:)
        integer(I8), allocatable, intent(out)   :: s1(:), s2(:)
        real(R8), allocatable, intent(out), optional    :: s1r(:), s2r(:)

        ! Auxiliary
        integer(I8)                             :: ni, counter, sz, &
            szold, szmult 
        real(R8)                                :: xe1, ye1, xe2, ye2, &
            d
        real(R8), allocatable                   :: tempx(:), tempy(:), &
            xi(:), yi(:), temps1r(:), temps2r(:)
        integer(I8), allocatable                :: temps1(:), &
            temps2(:), si(:), sortindex(:)
        logical, allocatable                    :: keepind(:)

        ! Loop
        integer(I8)                             :: i 

        ! Memory mgmt
        integer(I8), allocatable                :: mgmti(:, :)
        real(R8), allocatable                   :: mgmtr(:, :)

        ! Initialize
        !===========
        ! Checks
        if (allocated(x)) then
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if
        if (allocated(s1)) then 
            deallocate(s1) 
        end if
        if (allocated(s2)) then 
            deallocate(s2) 
        end if
        if (present(s1r)) then 
            if (allocated(s1r)) then 
                deallocate(s1r)
            end if
        end if 
        if (present(s2r)) then 
            if (allocated(s2r)) then 
                deallocate(s2r)
            end if 
        end if 

        ! Check if either both s1r, s2r are present or not present
        if (((.not. present(s1r)) .and. present(s1r)) &
            .or. (present(s1r) .and. .not. present(s2r))) then 
            call gdErrorHandler('PolygonIntersections: s1r and s2r should either be ' // &
                'both present or not present, one of the two is not supported')
        end if 

        ! Initialize
        counter     = 0 ! intersection counter 
        ni          = 0
        szold       = 0
        sz          = 2 ! initial size of intersection array
        szmult      = 2 ! size multiplier 

        ! Allocate
        allocate(tempx(sz), tempy(sz), temps1(sz), temps2(sz), &
            temps1r(sz), temps2r(sz))

        ! Compute intersections
        !======================
        ! Loop over p2
        !$omp parallel do default(none) schedule(dynamic) if(.not. omp_in_parallel()) &
        !$omp shared(x1, y1, x2, y2, counter, sz, tempx, tempy, temps1, &
        !$omp temps2, temps1r, temps2r, szmult) &
        !$omp private(i, xe1, ye1, xe2, ye2, xi, yi, si, ni, szold, &
        !$omp mgmti, mgmtr)
        do i = 1, size(x2)-1 
            ! Get coordinates of next polygon edge
            xe1 = x2(i)
            ye1 = y2(i)
            xe2 = x2(i+1)
            ye2 = y2(i+1)

            ! Compute intersections (duplicates of p1 with segment are 
            ! already removed in this routine)
            call SegmentSimplePolygonIntersections(x1, y1, xe1, ye1, xe2, ye2, &
                xi, yi, si)

            ! Check if intersections were found
            ni = size(xi)
            if (ni > 0) then 
                ! Memory MGMT
                !$omp critical
                if (counter + ni > sz) then 
                    
                    ! Store old size
                    szold = sz

                    ! Store old values
                    allocate(mgmti(szold, 2), mgmtr(szold, 2))
                    mgmti(:, 1) = temps1 
                    mgmti(:, 2) = temps2   
                    mgmtr(:, 1) = tempx 
                    mgmtr(:, 2) = tempy

                    ! Adjust size
                    do while (sz < counter+ni)
                        sz = sz*szmult 
                    end do

                    ! Reallocate
                    deallocate(tempx, tempy, temps1, temps2, temps1r, temps2r) 
                    allocate(tempx(sz), tempy(sz), temps1(sz), temps2(sz), &
                        temps1r(sz), temps2r(sz))

                    ! Add
                    tempx(1:szold) = mgmtr(:, 1) 
                    tempy(1:szold) = mgmtr(:, 2)
                    temps1(1:szold) = mgmti(:, 1)
                    temps2(1:szold) = mgmti(:, 2)

                    ! Deallocate mgmt arrays
                    deallocate(mgmti, mgmtr)
                end if 

                ! Add intersections
                tempx(counter+1:counter+ni) = xi 
                tempy(counter+1:counter+ni) = yi
                temps1(counter+1:counter+ni) = si 
                temps2(counter+1:counter+ni) = i

                ! Update counter
                counter = counter + ni
                !$omp end critical
            end if 
        end do  
        !$omp end parallel do

        ! Add to output
        allocate(x(counter), y(counter), s1(counter), s2(counter)) 
        x = tempx(1:counter) 
        y = tempy(1:counter) 
        s1 = temps1(1:counter) 
        s2 = temps2(1:counter) 

        ! Hedge for duplicates
        !=====================
        ! We only need to check s2 since duplicates of s1 should have
        ! already been removed before

        ! Sort intersections along s2 - may not be sorted anymore due to 
        ! parallellization
        allocate(sortindex(counter))
        call Sort(s2, ind=sortindex)
        x = x(sortindex)
        y = y(sortindex)
        s1 = s1(sortindex)

        ! Check which ones to keep
        allocate(keepind(counter))
        keepind = .true.
        do i = 1, counter-1
            ! Check for segment index 
            if ((s2(i+1) - s2(i)) == 1) then 
                ! Check for same intersection
                call Distance(d, x(i), y(i), x(i+1), y(i+1))
                if (d <= disttol) then 
                    ! Delete
                    keepind(i+1) = .false. 
                end if 
            end if 
        end do 

        ! Delete
        x = pack(x, keepind)
        y = pack(y, keepind)
        s1 = pack(s1, keepind)
        s2 = pack(s2, keepind)
        counter = count(keepind)

        ! Compute true intersection locations
        !====================================
        if (present(s1r) .and. present(s2r)) then 
            ! Compute the continuous intersection index (0: first point
            ! of polygon, ne+1: last point of polygon) - note: we need
            ! to use the vert array here instead of the edges array, since
            ! the latter is not necessarily sorted!
            allocate(s1r(counter), s2r(counter))

            ! First polygon index
            !s1r = ComputeI(x, y, p1%x(p1%edges(s1, 1)), p1%y(p1%edges(s1, 1)), &
            !    p1%x(p1%edges(s1, 2)), p1%y(p1%edges(s1, 2))) + s1 - 1 
            s1r = ComputeI(x, y, x1(s1), y1(s1), x1(s1+1), y1(s1+1)) + s1 - 1 

            ! Second polygon index
            !s2r = ComputeI(x, y, p2%x(p2%edges(s2, 1)), p2%y(p2%edges(s2, 1)), &
            !    p2%x(p2%edges(s2, 2)), p2%y(p2%edges(s2, 2))) + s2 - 1
            s2r = ComputeI(x, y, x2(s2), y2(s2), x2(s2+1), y2(s2+1)) + s2 - 1

        end if 

        ! Housekeeping
        !=============

        deallocate(tempx, tempy, temps1, temps2)


    end subroutine

    ! Inpolygon routine
    subroutine Inpolygon(polygon, xq, yq, in)

        ! Description
        !============
        ! Compute whether the query points xq, yq are located in the
        ! polygon or not. To compute this, we apply a winding number 
        ! algorithm (often referred to as Sunday's algorithm, but here 
        ! we based ourselves on the paper "A winding number and 
        ! point-in-polygon algorithm" by David G. Alciatore and Rick 
        ! Miranda. Basically, this consists of a ray tracing algorithm,
        ! where a ray is traced along the x-axis from the origin (the 
        ! polygon is shifted by each query point), and it is checked
        ! whether a polygon edge crosses the x-axis at the positive 
        ! side. If it goes from positive to negative, the winding number
        ! decreases by a half, if it goes from negative to positive y, 
        ! the winding number increases by a half. Points that lie 
        ! in the interior of the polygon should have a winding number
        ! different than zero. If they lay outside, it is equal to zero.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT), intent(in)       :: polygon 
        real(R8), allocatable, intent(in)   :: xq(:), yq(:) 
        logical, allocatable, intent(out)   :: in(:)

        ! Auxiliary 
        integer(I8)                         :: nq, np 
        integer(I8), allocatable            :: w(:), pvert(:) 

        real(R8)                            :: xi
        real(R8), allocatable               :: xp(:), yp(:), xcross(:), &
            tempx(:), tempy(:) 

        ! Loop
        integer(I8)                         :: i, j 

        ! Initialize
        !===========
        ! Associate
        associate(&
            x           => polygon%x,           &
            y           => polygon%y,           &
            vert        => polygon%vert,        &
            isclosed    => polygon%isclosed,    &
            ne          => polygon%ne) 

        ! Number of query points
        nq = size(xq)

        ! Winding number*2 (to be able to store as integer)
        allocate(w(nq), in(nq))
        w(:) =  0
        in(:) = .true.

        ! Shifted polygon points
        if (isclosed) then 
            allocate(xp(ne+1), yp(ne+1), pvert(ne+1))
            pvert = vert
        else
            allocate(xp(ne+2), yp(ne+2), pvert(ne+2))
            pvert = [vert, vert(1)]
        end if
        np = size(xp)
        allocate(xcross(np-1), tempx(np), tempy(np))
        tempx = x(pvert)
        tempy = y(pvert)
        
        ! Loop
        !=====
        do i = 1, nq 
            ! Shift points
            xp = tempx - xq(i)
            yp = tempy - yq(i)

            ! Compute which points cross x-axis
            xcross = yp(1:np-1)*yp(2:np)

            do j = 1, np-1
                if (xcross(j) < 0) then 
                    ! x-coordinate of intersection
                    xi = xp(j) + yp(j)*(xp(j+1) - xp(j))/(yp(j) - yp(j+1))

                    if (xi > 0) then 
                        if (yp(j) <0) then 
                            w(i) = w(i) + 2
                        else 
                            w(i) = w(i) - 2 
                        end if 
                    end if 
                elseif ( (yp(j) == 0) .and. (xp(j) > 0) ) then
                    if ( (yp(j+1) > 0) ) then 
                        w(i) = w(i) + 1
                    elseif (yp(j+1) == 0) then 
                        ! Do nothing - we're still on the same line and 
                        ! didn't cross
                    else 
                        w(i) = w(i) - 1
                    end if 
                elseif ( (yp(j+1) == 0) .and. (xp(j+1) > 0) ) then 
                    if (yp(j) < 0) then 
                        w(i) = w(i) + 1
                    elseif (yp(j) == 0) then 
                        ! Do nothing
                    else 
                        w(i) = w(i) - 1
                    end if 
                end if 

            end do

        end do 

        ! Check in/out
        !=============
        where (w == 0) in = .false.

        ! Housekeeping
        !=============
        deallocate(xp, yp, w, xcross)
        end associate 



    end subroutine 

    ! Edge overlap checker
    function CheckEdgeOverlap(x11, y11, x12, y12, x21, y21, &
        x22, y22) result(isoverlapping)

        ! Description
        !============
        ! This function checks whether two edges 'overlap', in the 
        ! sense that the boxes formed around these edges overlap. The 
        ! boxes have edges parallel with the axes. 

        ! Algorithm
        !==========
        ! For two boxes to not overlap, the common set of points of the 
        ! x-interval of both boxes should be the empty set (or the same
        ! for the y-interval). We hedge for distance precision tolerance
        ! as defined by macheps (i.e. we make the intervals disttol 
        ! larger on each side

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)        :: x11, y11, x12, y12, x21, y21, &
            x22, y22
        logical                     :: isoverlapping

        ! Auxiliary

        ! Loop

        ! Check boxes
        !============
        ! Initialize
        isoverlapping = .true. 

        ! x-interval
        if ((x11+disttol < x21-disttol) .and. (x11+disttol < x22-disttol) &
            .and. (x12+disttol < x21-disttol) .and. (x12+disttol < x22-disttol)) then 
                isoverlapping = .false. 
        elseif ((x21+disttol < x11-disttol) .and. (x21+disttol < x12-disttol) &
            .and. (x22+disttol < x11-disttol) .and. (x22+disttol < x12-disttol)) then 
                isoverlapping = .false. 
        elseif  ((y11+disttol < y21-disttol) .and. (y11+disttol < y22-disttol) &
            .and. (y12+disttol < y21-disttol) .and. (y12+disttol < y22-disttol)) then 
                isoverlapping = .false. 
        elseif ((y21+disttol < y11-disttol) .and. (y21+disttol < y12-disttol) &
            .and. (y22+disttol < y11-disttol) .and. (y22+disttol < y12-disttol)) then 
                isoverlapping = .false.    
        end if

       ! if ( (max(x11, x12)+disttol < min(x21, x22)-disttol) .or. &
       !     (max(x21, x22)+disttol < min(x11, x12)-disttol) ) then 
       !     isoverlapping = .false. 
       ! elseif ( (max(y11, y12)+disttol < min(y21, y22)-disttol) .or. &
       !     (max(y21, y22)+disttol < min(y11, y12)-disttol) ) then 
       !     isoverlapping = .false.
       ! end if 

    end function

    ! Edge overlap checker
    function CheckEdgeOverlap1D(x11, y11, x12, y12, x21, y21, &
        x22, y22) result(isoverlapping)

        ! Description
        !============
        ! This function checks whether two edges 'overlap', in the 
        ! sense that the boxes formed around these edges overlap. The 
        ! boxes have edges parallel with the axes. 

        ! Algorithm
        !==========
        ! For two boxes to not overlap, the common set of points of the 
        ! x-interval of both boxes should be the empty set (or the same
        ! for the y-interval). We hedge for distance precision tolerance
        ! as defined by macheps (i.e. we make the intervals disttol 
        ! larger on each side

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)  :: x11, y11, x12, y12
        real(R8), intent(in)                :: x21, y21, x22, y22
        logical, allocatable, dimension(:)  :: isoverlapping

        ! Auxiliary

        ! Loop
        integer(I8)                         :: i 

        ! Check boxes
        !============
        ! Initialize
        allocate(isoverlapping(size(x11)))
        isoverlapping = .true. 

        ! x-interval
        do i = 1, size(x11)
                if ((x11(i)+disttol < x21-disttol) .and. (x11(i)+disttol < x22-disttol) &
                .and. (x12(i)+disttol < x21-disttol) .and. (x12(i)+disttol < x22-disttol)) then 
                    isoverlapping(i) = .false. 
            elseif ((x21+disttol < x11(i)-disttol) .and. (x21+disttol < x12(i)-disttol) &
                .and. (x22+disttol < x11(i)-disttol) .and. (x22+disttol < x12(i)-disttol)) then 
                    isoverlapping(i) = .false. 
            elseif  ((y11(i)+disttol < y21-disttol) .and. (y11(i)+disttol < y22-disttol) &
                .and. (y12(i)+disttol < y21-disttol) .and. (y12(i)+disttol < y22-disttol)) then 
                    isoverlapping(i) = .false. 
            elseif ((y21+disttol < y11(i)-disttol) .and. (y21+disttol < y12(i)-disttol) &
                .and. (y22+disttol < y11(i)-disttol) .and. (y22+disttol < y12(i)-disttol)) then 
                    isoverlapping(i) = .false.    
            end if
        end do 
        
        !where ((x11+disttol < x21-disttol) .and. (x11+disttol < x22-disttol) &
        !    .and. (x12+disttol < x21-disttol) .and. (x12+disttol < x22-disttol))  
        !        isoverlapping = .false. 
        !end where 
        !where ((x21+disttol < x11-disttol) .and. (x21+disttol < x12-disttol) &
        !    .and. (x22+disttol < x11-disttol) .and. (x22+disttol < x12-disttol))  
        !        isoverlapping = .false. 
        !end where 
        !where  ((y11+disttol < y21-disttol) .and. (y11+disttol < y22-disttol) &
        !    .and. (y12+disttol < y21-disttol) .and. (y12+disttol < y22-disttol))  
        !        isoverlapping = .false. 
        !end where 
        !where ((y21+disttol < y11-disttol) .and. (y21+disttol < y12-disttol) &
        !    .and. (y22+disttol < y11-disttol) .and. (y22+disttol < y12-disttol))  
        !        isoverlapping = .false.    
        !end where

       ! if ( (max(x11, x12)+disttol < min(x21, x22)-disttol) .or. &
       !     (max(x21, x22)+disttol < min(x11, x12)-disttol) ) then 
       !     isoverlapping = .false. 
       ! elseif ( (max(y11, y12)+disttol < min(y21, y22)-disttol) .or. &
       !     (max(y21, y22)+disttol < min(y11, y12)-disttol) ) then 
       !     isoverlapping = .false.
       ! end if 

    end function

    ! Edge coincidence checker
    subroutine CheckSegmentCoincidence(x, y, x11, y11, x12, y12, x21, &
        y21, x22, y22)

        ! Description
        !============
        ! Check if two segments are coincident, either at one of the two
        ! end vertices or at both. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)    :: x11, y11, x12, y12, x21, y21, &
            x22, y22 
        real(R8), intent(out)   :: x, y

        ! Auxiliary
        real(R8)                :: nan, d11, d12, d21, d22


        ! Compute distances
        !==================
        ! Distances between segment points
        call Distance(d11, x11, y11, x21, y21)
        call Distance(d12, x11, y11, x22, y22)
        call Distance(d21, x12, y12, x21, y21)
        call Distance(d22, x12, y12, x22, y22)

        ! Check
        if (( ((d11 < disttol) .and. (d22 < disttol)) &
                .or. ((d12 < disttol) .and. (d21 < disttol)) )) then 
            ! Edges are the same up to distance tolerance, take
            ! average
            x = 0.5*(x11 + x12)
            y = 0.5*(y11 + y12)
        else
            
            ! Checks
            if (d11 < disttol) then 
                ! First points coincide
                x = x11;
                y = y11;
            elseif (d12 < disttol) then 
                ! First point of first segment and second point of
                ! second segment coincide
                
                x = x11;
                y = y11;
            elseif (d21 < disttol) then 
                ! Second point of first segment and first point of
                ! second segment coincide
                
                ! Intersection in single point
                x = x12;
                y = y12;
            elseif (d22 < disttol) then 
                ! Second point of first segment and first point of
                ! second segment coincide
                
                x = x12;
                y = y12;
            else
                ! No coincidence of nodes, set to NaN
                x = IEEE_VALUE(nan, IEEE_QUIET_NAN)
                y = IEEE_VALUE(nan, IEEE_QUIET_NAN)
            end if 
        end if


    end subroutine

    ! Edge inbetween edges checker
    function IsEdgeInClosedSimplePolygon(xp, yp, vp, v1, v2) &
        result(isinbetween)

        ! Description
        !============
        ! This function computes whether the vertex pair, defined by
        ! indices [v1, v2], lies in the interior of the simple 
        ! (closed) polygon. It is assumed that start and end point
        ! are the same vertex. Note that if the vertex pair consists of 
        ! consecutive indices, the edge lies automatically in the 
        ! polygon. If this is not the case, then we need to check 
        ! the following:
        !
        ! - does the edge start and end in the interior of the polygon?
        !   if no -> edge lies outside, if yes -> continue
        ! - does the edge intersect with any non-neighbouring edge? 
        !   if no -> edge lies in interior, if yes -> edge lies outside
        !
        ! Note that this routine naturally only supports edges formed
        ! by polygon vertices. The inputs should be the coordinates of 
        ! the (subsequent) edges, and vp the vertex ID. The latter 
        ! should be used to check if vertices are the same in order to 
        ! allow polygons that have 'touching' boundaries. 

        ! Note: it is assumed that the polygon is not self-intersecting.
        ! Otherwise, surface area computations may be off and the result
        ! will be wrong.

        ! Note: the first check on starting/ending in the interior of 
        ! the domain is done based on:
        ! - the angle between neighbouring edges of the first edge vertex, 
        !   measured in the interior of the polygon (alpha)
        ! - the angle between the first neighbouring edge (in vertex 
        !   order of the polygon) and the edge vector, pointing away
        !   from the first vertex (beta)
        ! - similar angles, but for the second vertex and its 
        !   neighbouring edges
        ! These angles should always be in the interval [0, 2*pi]. The
        ! edge then lies within the polygon if for both vertices it 
        ! holds that:
        !
        !       0 <= beta <= alpha

        ! Note: v1 and v2 are indices in xp, yp, vp

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)          :: xp, yp
        integer(I8), intent(in)                     :: v1, v2, vp(:)
        logical                                     :: isinbetween 

        ! Auxiliary
        integer(I8)                                 :: nv, v1s, v1e, &
            v2s, v2e
        real(R8)                                    :: dx, dy, &
            dxs1, dys1, dxe1, dye1, dxs2, dys2, dxe2, dye2, & 
            d0, xint, yint, alpha(1:2), beta(1:2), &
            dpalpha(1:2), cpalpha(1:2), dpbeta(1:2), cpbeta(1:2), &
            surfA

        ! Loop
        integer(I8)                                 :: i 

        ! Initialize
        !===========
        ! Output
        isinbetween = .false. 

        ! Check dimensions
        if ((size(xp) /= size(yp)) .or. (size(xp) /= size(vp))) then 
            call gdErrorHandler('IsEdgeInClosedSimplePolygon: inputs ' // & 
                'have inconsistent dimensions')
        end if 
        nv = size(xp)

        ! Check last point
        call Distance(d0, xp(1), yp(1), xp(nv), yp(nv))
        if ((d0 > disttol) .and. (vp(1) /= vp(nv))) then 
            call gdErrorHandler('IsEdgeInClosedSimplePolygon: first and ' // &
                'last vertex should be the same, but have non-zero distance')
        end if 

        ! Compute
        !========
        ! Check for trivial case
        if ((v1 == v2) .or. (abs(v1 - v2) == 1) .or. &
            (v1 == 1 .and. v2 == nv-1) .or. (v1 == nv-1 .and. v2 == 1)) then 
            ! All cases where the edge lies on the polygon
            isinbetween = .true.
            return 
        end if 

        ! Get vertices of previous and next vertex
        if (v1 /= 1) then 
            v1s = v1 - 1
        else
            v1s = nv-1
        end if 
        if (v1 /= nv) then 
            v1e = v1 + 1
        else
            v1e = 2
        end if 
        if (v2 /= 1) then 
            v2s = v2 - 1
        else
            v2s = nv-1
        end if 
        if (v2 /= nv) then 
            v2e = v2 + 1
        else
            v2e = 2
        end if 

        ! Unpack
        associate(&
            x1      => xp(v1),  &
            y1      => yp(v1),  &
            x1s     => xp(v1s), &
            y1s     => yp(v1s), &
            x1e     => xp(v1e), &
            y1e     => yp(v1e), &
            x2s     => xp(v2s), &
            y2s     => yp(v2s), &
            x2e     => xp(v2e), &
            y2e     => yp(v2e), &
            x2      => xp(v2),  &
            y2      => yp(v2)   &
            )

        ! Compute distances along the polygon direction (so from vert 1
        ! to vert N)
        dx = x2 - x1 
        dy = y2 - y1
        dxs1 = (x1s - x1 )
        dys1 = (y1s - y1)
        dxe1 = x1e - x1 
        dye1 = y1e - y1
        dxs2 = (x2s - x2) 
        dys2 = (y2s - y2)
        dxe2 = x2e - x2 
        dye2 = y2e - y2

        ! Compute dot and cross products
        dpalpha(1) = (dxs1*dxe1 + dys1*dye1)
        dpalpha(2) = (dxs2*dxe2 + dys2*dye2)
        cpalpha(1) = dxs1*dye1 - dys1*dxe1
        cpalpha(2) = dxs2*dye2 - dys2*dxe2
        dpbeta(1) = (dxs1*dx + dys1*dy)
        dpbeta(2) = -(dxs2*dx + dys2*dy)
        cpbeta(1) = (dxs1*dy - dys1*dx)
        cpbeta(2) = -(dxs2*dy - dys2*dx)

        ! Compute angles
        alpha = atan2(cpalpha, dpalpha)
        beta = atan2(cpbeta, dpbeta)
        !where (alpha < 0.0_R8) alpha = alpha + 2*pi_R8
        !where (beta < 0.0_R8) beta = beta + 2*pi_R8

        ! Compute surface area
        surfA = ComputeSimplePolygonSurfaceArea(xp, yp)

        ! Adjust angles
        if (surfA >= 0.0_R8) then 
            ! If angle is negative, then concave angle
            where (alpha < 0.0_R8) alpha = alpha + 2*pi_R8 
            where (beta < 0.0_R8)  beta = beta + 2*pi_R8
            if (all(beta <= alpha)) then 
                isinbetween = .true.
            else
                isinbetween = .false.
            end if 
        else
            ! If angle is negative, then convex angle
            where (alpha < 0.0_R8) alpha = alpha + 2*pi_R8 
            where (beta < 0.0_R8)  beta = beta + 2*pi_R8
            if (all(beta >= alpha)) then 
                isinbetween = .true.
            else
                isinbetween = .false.
            end if 
        end if 

        ! First check
        if (.not. isinbetween) then 
            return 
        end if 

        ! Check intersections with non-neighbouring edges
        do i = 1, nv-1 
            ! Next edge is [i, i+1]
            ! Check if we should skip
            if ((vp(i) == vp(v1)) .or. (vp(i+1) == vp(v1)) .or. &
                (vp(i) == vp(v2)) .or. (vp(i+1) == vp(v2)) .or. &
                ((v1 == 1) .and. (i == nv .or. i+1 == nv)) .or. & 
                ((v2 == 1) .and. (i == nv .or. i+1 == nv)) &
                ) then 
                cycle 
            end if 

            ! Compute intersections
            call SegmentIntersections(xint, yint, xp(i), yp(i), &
                xp(i+1), yp(i+1), x1, y1, x2, y2)

            ! Check
            if (.not. isnan(xint)) then
                isinbetween = .false. 
                exit 
            end if 

        end do

        ! Housekeeping
        end associate

    end function

    ! Polygon surface area
    function ComputePolygonSurfaceArea(polygon) result(out)

        ! Description
        !============
        ! Simple function that computes the (signed) surface area 
        ! enclosed by the polygon. For this, we employ the simple
        ! trapezoidal rule for integration, which is exact for 
        ! piecewise-linear polygons. The surface area is negative
        ! in case of a right-turning polygon, and positive for a
        ! left-turning.

        ! Note: this is the closed polygon surface area. If the polygon
        ! is open, the last vertex will be connected to the first one
        ! for the area computation.

        ! Declare variables
        !==================
        ! Arguments
        class(PolygonUDT)                       :: polygon 
        real(R8)                                :: out 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: dx, yf 

        ! Compute
        !========
        dx = polygon%x(polygon%vert) - polygon%x([polygon%vert(2:polygon%nv), polygon%vert(1)])
        yf = polygon%y(polygon%vert) + polygon%y([polygon%vert(2:polygon%nv), polygon%vert(1)])
        out = 0.5*sum(-dx*yf)



    end function

    function ComputeSimplePolygonSurfaceArea(x, y) result(out)

        ! Description
        !============
        ! Same as function for polygon, but now with simple x, y input

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in), dimension(:)      :: x, y
        real(R8)                                :: out 

        ! Auxiliary
        real(R8), allocatable, dimension(:)     :: dx, yf 

        ! Compute
        !========
        dx = x - [x(2:size(x)), x(1)]
        yf = y + [y(2:size(y)), y(1)]
        out = 0.5*sum(-dx*yf)

    end function

    !------------------------------------------------------------------!
    !                               Writing                            !
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    !                                Tests                             !
    !------------------------------------------------------------------!

    ! Main test routine
    subroutine TestPolygonRoutines()

        ! Description
        !============
        ! Run tests on some polygons

        ! Declare variables
        !==================
        ! Arguments

        ! Auxiliary
        type(PolygonSetUDT)         :: ps 

        integer(I8)                 :: id 
        real(R8), allocatable       :: x(:), y(:) 


        ! Set polygon
        !============
        ! Polygon case
        id = 1

        ! Coordinates
        call ConstructTestPolygonCoordinates(id, x, y)

        ! Test
        !=====
        call ps%Construct(x, y)

        call Plot2DPolygon(ps%polygons(1)%x, ps%polygons(1)%y, size(ps%polygons(1)%x), '-p')
        call Quiverplot2D(ps%polygons(1)%ex, ps%polygons(1)%ey, &
            ps%polygons(1)%nx, ps%polygons(1)%ny, &
            size(ps%polygons(1)%ex), '-p')

        ! Set polygon
        !============
        ! Polygon case
        id = 5

        ! Coordinates
        call ConstructTestPolygonCoordinates(id, x, y)

        ! Test
        !=====
        if (allocated(ps%polygons)) then 
            deallocate(ps%polygons)
        end if
        call ps%Construct(x, y)

    end subroutine

    ! Set of test polygons
    subroutine ConstructTestPolygonCoordinates(id, x, y) 

        ! Description
        !============
        ! Return the test polygon coordinates belonging to index id

        ! Declare variables
        !==================
        ! Arguments
        integer(I8), intent(in)         :: id 
        real(R8), allocatable           :: x(:), y(:)

        ! Auxiliary
        integer(I8)                     :: np

        ! Construct polygon coordinates
        !==============================
        ! Check
        if (allocated(x)) then 
            deallocate(x) 
        end if 
        if (allocated(y)) then 
            deallocate(y) 
        end if 

        ! Build
        select case (id) 

            case (1) 

                ! Normal polygon
                np = 4
                allocate(x(np), y(np))
                x = [1, 1, 2, 4]
                y = [0, 2, 2, 1]
            
            case (5)

                ! Polygon with some duplicate points
                np = 12
                allocate(x(np), y(np))
                x = [1, 1, 1, 1, 2, 4, 1, 1, 2, 5, 6, 8]
                y = [1, 1, 0, 2, 2, 1, 1, 1, 3, 2, 0, 2]

            case default 

            ! Throw error
            call PolygonErrorHandler('No coordinates found for this id')
            
        end select 




    end subroutine



end module