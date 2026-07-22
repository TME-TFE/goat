!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains utility routines for the optimization modules of
! the grid deformation framework. These routines include for example 
! routines to compute faces on which orthogonality constraints should be
! imposed, routines that allow to derive the vertices that are x-points,
! and so on. In general, any routine that is necessary to set up the 
! constraints/costfunction, but that is not needed anymore for the 
! remainder of the code, can be added here. 

module gdmod_utility_optimization

    ! Initialize
    !============
    ! Load modules
    !use mod_plotter
    use gdmod_types
    use gdmod_userinput 
    use mod_definitions, only: targetID, vesselID, outerboundaryID, &
        coreID, interiorID
    !use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    ! The usual
    implicit none
    save
    public 

    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Grid quantities                       !
    !------------------------------------------------------------------!

    ! Determination of X-point(s)
    subroutine DetermineXPoints(xpind, nxpind, order, grid)

        ! Description
        !============
        ! This routine determines the indices (and their number) based 
        ! on the information given in the grid structure. Note that 
        ! the location of the x-point should already be available 
        ! through output of the grid generator, but we provide a routine 
        ! here to recompute the location if it wasn't saved explicitly
        ! in the grid structure. Additionally, we determine the order of
        ! the x-point (see algorithm description below)

        ! Algorithm
        !==========
        ! Basically, we simply check if a vertex has 4 or more 
        ! neighbours with the same flux surface ID (which is non-zero)
        ! as the current vertex. The order is then simply determined by 
        ! the number of separatrix legs, i.e. 
        !
        !       o = n/2 - 1,
        !
        ! where n is the number of vertices with the same ID and o is 
        ! the order. 

        ! Initialize
        !===========
        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8)                         :: nxpind
        integer(I8), allocatable            :: xpind(:), order(:)
        type(GridUDT), intent(in)           :: grid

        ! Loop variables
        integer(I8)                         :: i

        ! Auxiliary variables 
        integer(I8)                         :: tfID, ncIDs
        integer(I8), allocatable            :: temporder(:), &
            tempxpind(:), tvn(:), tvnfID(:)

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            fID         => grid%vert%fieldlineID,   &
            vert        => grid%vert)

        ! Check allocation - shouldn't be the case as size unknown a 
        ! priori
        if (allocated(xpind)) then
            ! Deallocate
            deallocate(xpind)
        end if
        if (allocated(order)) then 
            ! Deallocate
            deallocate(order)
        end if

        ! Allocate temporary arrays (too big, trim later)
        nxpind = 0
        allocate(tempxpind(vert%ntot))
        allocate(temporder(vert%ntot))

        ! Determine x-points
        !===================
        ! Loop over all vertices
        do i = 1, vert%ntot 
            ! Get the current field line ID, skip if zero
            tfID = fID(i)
            if (tfID == 0) then
                cycle ! skip rest of the loop for this index
            end if

            ! Get the vertex neighbours
            allocate(tvn(vert%neigP(i, 2)))
            tvn = vert%neig(vert%neigP(i, 1):&
                (vert%neigP(i, 1) + vert%neigP(i, 2)-1))

            ! Get their IDs
            allocate(tvnfID(vert%neigP(i, 2)))
            tvnfID = fID(tvn)

            ! Get the number of common IDs
            ncIDs = count(tvnfID == tfID)

            ! If equal or larger than 4, add as x-point
            if (ncIDs >= 4) then 
                ! Check if it is a multiple of 2
                if (modulo(ncIDs, 2) .ne. 0) then 
                    ! Uneven number, throw error
                    call gdErrorHandler('DetermineXPoints: '&
                    // 'supposed x-point has uneven number' &
                    // 'of vertices with same ID, check grid' &
                    // ' consistency')
                end if

                ! Update the counter
                nxpind = nxpind + 1

                ! Add the x-point
                tempxpind(nxpind) = i 

                ! Compute the order
                temporder(nxpind) = ncIDs/2 - 1

            end if
            
            ! Deallocate
            deallocate(tvn, tvnfID)

        end do

        ! Output
        !=======
        ! Allocate
        allocate(xpind(nxpind), order(nxpind))

        ! Set output
        xpind = tempxpind(1:nxpind)
        order = temporder(1:nxpind)

        ! Housekeeping
        !=============
        ! Deallocate 
        deallocate(tempxpind, temporder)

        ! Deassociate
        end associate


    end subroutine

    ! Determination of X-point edges that are aligned with flux 
    ! surfaces
    subroutine DetermineFluxAlignedXPointEdges(nxpointedges, &
        xpointedges, grid)

        ! Description
        !============
        ! This routine returns the edges (as vertex pairs) that have an
        ! x-point and a neighbour of that x-point with the same flux 
        ! surface ID. 

        ! Algorithm
        !==========
        ! 1) Determine the x-points by calling the DetermineXPoints 
        !   routine
        ! 2) For each x-point, determine which of its vertex neighbours
        !   has the same flux ID as the x-point. 

        ! Notes
        !======
        ! Note 1: since the order and number of the xpoints are known, 
        ! the total number of edges can be determined a priori. 

        ! Initialize

        ! Initialize
        !===========
        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8)                         :: nxpointedges
        integer(I8), allocatable            :: xpointedges(:, :)
        type(GridUDT), intent(in)           :: grid

        ! Loop variables
        integer(I8)                         :: i, ec 

        ! Auxiliary variables 
        integer(I8)                         :: ntvn, nxpind

        integer(I8), allocatable            :: xpind(:), order(:), &
            tvn(:)

        logical, allocatable                :: isEqualID(:)

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => grid%vert)

        ! Get x-points
        call DetermineXPoints(xpind, nxpind, order, grid)

        ! Compute number of edges
        nxpointedges = sum(order*2 + 2)
        
        ! Allocate (soft)
        if (allocated(xpointedges)) then 
            print *, 'DetermineFluxAligneXPointEdges: ' &
                // 'reallocating xpointedges'
            deallocate(xpointedges)
        end if
        allocate(xpointedges(nxpointedges, 2))

        ! Get x-point edges
        !==================
        ! Initialize edge counter
        ec = 0

        ! Loop over x-points
        do i = 1, nxpind 
            ! Get the number of neighbours
            ntvn = vert%neigP(xpind(i), 2)

            ! Get neighbours
            allocate(tvn(ntvn))
            tvn = vert%neig(vert%neigP(xpind(i), 1):&
                vert%neigP(xpind(i), 1)+ntvn-1)
            
            ! Check which vertices have the same flux surface ID
            allocate(isEqualID(ntvn))
            isEqualID = vert%fieldlineID(xpind(i)) == &
                vert%fieldlineID(tvn)
            
            ! Sanity check: should be equal amount as order*2 + 2
            if (count(isEqualID) .ne. (2*order(i)+2) ) then 
                ! Print vertex ID of current x-point
                print *, 'X-point vertex ID: ', xpind(i)

                ! Throw error
                call gdErrorHandler('DetermineFluxAlignedXPointEdges:' &
                    // 'inconsistent amount of edges detected for ' &
                    // 'x-point with above mentioned vertex ID')
                
            end if

            ! Add the edges
            xpointedges(ec+1:ec+2*order(i)+2, 1) = xpind(i)
            xpointedges(ec+1:ec+2*order(i)+2, 2) = pack(tvn, isEqualID)

            ! Housekeeping
            deallocate(tvn)
        end do
        

        ! End associate
         end associate 


    end subroutine

    ! Determination of tangency points
    subroutine DetermineTangencyPoints(tpind, ntpind, tptype, grid)

        ! Description
        !============
        ! Determine the tangency points in the grid based on the grid
        ! topology. A tangency point is assumed to adhere to the 
        ! following criteria:
        ! 
        !   - Vertex lies on boundary
        !   - Vertex lies on a flux surface
        !   - Vertex has two vertex neighbours that have the same 
        !   flux surface ID
        !
        !   OR
        !
        !   - Vertex is a boundary vertex with unique flux surface ID
        !   - Vertex has only neighbours with non-zero flux surface ID
        !   that is the same for all vertices
        !
        ! Note that both definitions are required to capture all 
        ! tangency points, but that possibly only the first type
        ! of tangency points is important.

        ! Declare
        !========
        ! Arguments
        integer(I8)                     :: ntpind
        integer(I8), allocatable        :: tpind(:), tptype(:)
        type(GridUDT)                   :: grid 

        ! Auxiliary 
        integer(I8), allocatable        :: tv(:), tvneig(:), &
            tvneigID(:)
        logical, allocatable            :: istpID(:), isvesselvertex(:), &
            isvesselface(:), isuniqueID(:), tvec(:)

        ! Loop
        integer(I8)                     :: i 


        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => grid%vert,   &
            ID          => grid%vert%fieldlineID & 
            )

        ! Initialize
        allocate(istpID(grid%vert%ntot))
        istpID(:) = .false.

        ! Check
        if (allocated(tpind)) then 
            deallocate(tpind)
        end if 
        if (allocated(tptype)) then 
            deallocate(tptype)
        end if 
        
        ! Allocate
        allocate(tptype(grid%vert%ntot))
        tptype(:) = 0

        ! First type
        !===========
        ! Get vertices that lie on a flux surface and that are vessel
        ! vertices
        call DetermineVesselVertices(isvesselvertex, isvesselface, grid)
        tv = pack([(i, i=1, vert%ntot)], (isvesselvertex .and. (ID .ne. 0)))

        ! Check vertex neighbours
        do i = 1, size(tv, 1)
            ! Get neighbours
            tvneig = GetVertNeig(vert, tv(i))

            ! Check
            tvneigID = ID(tvneig)
            if (count(tvneigID == ID(tv(i))) == 2) then 
                istpID(tv(i)) = .true.
                tptype(tv(i)) = 1
            end if 
        end do

        ! Second type
        !============
        ! Initialize
        allocate(isuniqueID(vert%ntot), tvec(vert%ntot))
        isuniqueID(:) = .false.
        do i = 1, maxval(ID)
            tvec = grid%vert%fieldlineID == i
            if (count(tvec) == 1) then 
                where (tvec) isuniqueID = .true.
            end if 
        end do
        tv = pack([(i, i=1, vert%ntot)], isvesselvertex .and. (isuniqueID .or. grid%vert%fieldlineID == 0))

        ! Check vertex neighbours
        do i = 1, size(tv, 1)
            ! Get neighbours
            tvneig = GetVertNeig(vert, i)

            ! Check
            tvneigID = ID(tvneig)
            if (all(tvneigID == tvneigID(1) .and. (tvneigID(1) .ne. ID(tv(i))))) then 
                istpID(tv(i)) = .true.
                tptype(tv(i)) = 2
            end if  
        end do

        ! Set tpind
        !==========
        ntpind = count(istpID)
        allocate(tpind(ntpind))
        tpind = 0
        tpind = pack([(i, i = 1, vert%ntot)], istpID)
        tptype = pack(tptype, istpID)
    

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Determination of vessel edges that are aligned with flux surfaces
    subroutine DetermineFluxAlignedVesselEdges(nvesseledges, &
        vesseledges, grid, doTP, doWG)

        ! Description
        !============
        ! This routine returns the edges (as vertex pairs) that have one
        ! vessel node and a node that has the same flux surface ID as 
        ! the vessel node. Normally, this is one edge per edge that 
        ! lies on the vessel. This routine's intended use is for the 
        ! edge length constraints, where one may want to constrain the 
        ! cell size near the target, which is effectively achieved by 
        ! constraining the edge length. 
        !
        ! It is possible through the input argument logicals 'doTP' and 
        ! 'doWG' to include (or not) the target plates and general wide 
        ! grid (WG) boundaries as desired. 

        ! Algorithm
        !==========
        ! The algorithm is pretty straightforward. 
        ! First, we check which boundaries should be considered. For
        ! each boundary, we loop over the vertices that belong to that 
        ! boundary and determine which edge (if there is one) complies
        ! to the criteria mentioned above (i.e. same flux surface ID). 

        ! Notes
        !======
        ! Note 1: it is assumed that the target plate indices are either
        ! 1 or 2 (1 for inner, 2 for outer, though doesn't really matter
        ! here.)

        ! Note 2: we explicitly check that any neighbour of a vertex is
        ! not part of the current boundary as well. Otherwise, this 
        ! could give issues in some pathological cases (e.g. single 
        ! edge flux surfaces that start and end at the same boundary)

        ! Initialize
        !===========
        ! The usual
        implicit none 

        ! Declare variables
        !==================
        ! Arguments
        integer(I8)                         :: nvesseledges
        integer(I8), allocatable            :: vesseledges(:, :)
        logical, intent(in)                 :: doTP, doWG
        type(GridUDT), intent(in)           :: grid

        ! Loop variables
        integer(I8)                         :: i, j

        ! Auxiliary variables 
        integer(I8)                         :: ntv, ntvn, ne

        integer(I8), allocatable            :: tv(:), tvn(:), &
            tempvesseledges(:, :)

        logical                             :: condition

        logical, allocatable                :: mask(:), &
            isBndVertex(:), isEqualID(:), isVertIncluded(:)

        ! Data

        ! Initialize
        !===========
        ! Associate
        associate(&
            vert        => grid%vert,   &
            faces       => grid%face,  &
            bnd         => grid%bnd)

        ! Allocate (too big, trim later)
        allocate(tempvesseledges(faces%ntot, 2))
        allocate(isBndVertex(vert%ntot), isVertIncluded(vert%ntot))
            
        ! Initialize
        nvesseledges = 0
        isVertIncluded = .false. 

        ! Loop over all boundaries
        !=========================
        do i = 1, size(bnd) 
            ! Check if this boundary is a target plate or wide grid bnd
            condition = (bnd(i)%ID == targetID) .and. doTP
            condition = condition .or. ( (bnd(i)%ID == vesselID) .and. doWG)
            if (condition) then 
                ! Get the vertices of this boundary
                ntv = bnd(i)%nvert
                allocate(tv(ntv))
                tv = bnd(i)%vert 

                ! Set these vertices as boundary vertices in order to 
                ! exclude them in the edge determination
                isBndVertex(:) = .false.
                isBndVertex(tv) = .true.

                ! For each vertex, get the neighbours
                do j = 1, ntv
                    ! Skip if vertex was already included
                    if (isVertIncluded(tv(j))) then 
                        print *, 'vertex already included, skipping...'
                        cycle
                    end if 

                    ! Extract neighbours of this vertex
                    ntvn = vert%neigP(tv(j), 2)
                    allocate(tvn(ntvn))
                    tvn = vert%neig(vert%neigP(tv(j), 1):vert%neigP(tv(j), 1)+ntvn-1)

                    ! Initialize mask for inclusion
                    allocate(mask(ntvn))
                    mask(:) = .true.

                    ! Don't include any vertices that are now boundary
                    ! vertices
                    where (isBndVertex(tvn)) mask = .false. 

                    ! Don't include vertices that have a different field
                    ! line ID
                    allocate(isEqualID(ntvn))
                    isEqualID = (vert%fieldlineID(tv(j)) == &
                        vert%fieldlineID(tvn))
                    where (.not. isEqualID) mask = .false.

                    ! Add to temporary edges
                    ne = count(mask)
                    tempvesseledges(nvesseledges+1:nvesseledges+ne, 1) = &
                        tv(j)
                    tempvesseledges(nvesseledges+1:nvesseledges+ne, 2) = &
                        pack(tvn, mask)

                    ! Update counter
                    nvesseledges = nvesseledges + ne

                    ! Update logical
                    isVertIncluded(tv(j)) = .true.

                    ! Housekeeping
                    deallocate(tvn, mask, isEqualID)                    
                end do

                ! Deallocate
                deallocate(tv)

            end if
            
        end do

        ! Assign output
        !==============
        ! Soft allocation
        if (allocated(vesseledges)) then 
            print *, 'DetermineVesselEdges: reallocating vessel edges'
            deallocate(vesseledges)
        end if
        allocate(vesseledges(nvesseledges, 2))

        ! Assign
        vesseledges = tempvesseledges(1:nvesseledges, :)

        ! Housekeeping
        !=============
        ! Deallocate
        deallocate(isBndVertex, tempvesseledges)

        ! End associate
        end associate

    end subroutine

    ! Determine which vertices lie on the vessel 
    subroutine DetermineVesselVertices(isvesselvertex, isvesselface, grid)

        ! Description
        !============
        ! Determine which vertices lie on vessel boundaries. This is
        ! checked by the boundary identifiers. The following vessel
        ! boundaries are considered:
        ! - target plates (TP)
        ! - far vessel boundaries (V)
        ! The boolean 'isvesselvertex' and 'isvesselface' are returned

        ! Declare variables
        !==================
        ! Arguments
        logical, allocatable            :: isvesselvertex(:), &
            isvesselface(:)
        type(GridUDT)                   :: grid 

        ! Auxiliary
        integer(I8), allocatable        :: bffaces(:)

        ! Loop
        integer(I8)                     :: ib, i

        ! Initialize
        !===========
        ! Associate
        associate(&
            bnd     => grid%bnd & 
            )
        
        ! Initialize
        if (allocated(isvesselface)) then 
            deallocate(isvesselface)
        end if 
        if (allocated(isvesselvertex)) then 
            deallocate(isvesselvertex)
        end if 
        allocate(isvesselface(grid%face%ntot), &
            isvesselvertex(grid%vert%ntot))
        isvesselface(:)     = .false.
        isvesselvertex(:)   = .false.
        
        ! Vertices
        !=========
        do ib =  1, size(bnd)
            ! Check if boundary belongs to vessel boundaries
            select case (bnd(ib)%ID) 
                
            case (targetID, vesselID)

                isvesselvertex(bnd(ib)%vert) = .true.

            case default 

            end select
        end do 

        ! Faces
        !======
        ! Determine boundary faces
        allocate(bffaces(count(grid%face%BF)))
        bffaces = pack([(i, i=1,grid%face%ntot)], grid%face%BF)

        ! Check which boundary faces have two boundary vertices as
        ! determined above
        do i = 1, size(bffaces, 1)
            if (all(isvesselvertex(grid%face%vert(bffaces(i), :)))) then 
                isvesselface(bffaces(i)) = .true.
            end if 
        end do 
        
        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! Determine which flux surfaces are closed/open
    function DetermineFluxSurfaceClosure(grid) result(closuretype)

        ! Description
        !============
        ! This function returns all the flux surface IDs of flux 
        ! surfaces that are closed (1)/open(0)/undetermined(-1). 
        ! This is determined by checking if 
        ! all vertices occur exactly twice in the flux surface ID. 
        ! It is assumed that flux surface IDs go from 1 to 
        ! max(vert%fieldlineID) (flux surface IDs with no vertices 
        ! are not considered closed or open)

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(in)               :: grid 
        integer(I8), allocatable, dimension(:)  :: closuretype 

        ! Auxiliary
        integer(I8)                             :: maxfsID
        integer(I8), allocatable, dimension(:)  :: vertcount, &
            fsfaces

        ! Loop 
        integer(I8)                             :: i, j

        ! Initialize
        !===========
        ! Determine max amount of flux surfaces
        maxfsID = maxval(grid%vert%fieldlineID)

        ! Initialize output 
        allocate(closuretype(maxfsID))
        closuretype = 0_I8 

        ! Determine closure
        !==================
        ! Loop
        allocate(vertcount(grid%vert%ntot))
        do i = 1, maxfsID 
            ! Get the flux surface faces, if none -> undefined
            fsfaces = GetFSFace(grid%data%fluxdata, i)
            if (size(fsfaces) == 0) then 
                closuretype(i) = -1_I8 
                cycle 
            end if 

            ! Count vertex occurrences
            vertcount = 0_I8
            do j = 1, size(fsfaces)
                vertcount(grid%face%vert(fsfaces(j), :)) = &
                    vertcount(grid%face%vert(fsfaces(j), :)) + 1
            end do 

            ! Check
            if (any(vertcount == 1)) then 
                closuretype(i) = 0_I8
            else
                closuretype(i) = 1_I8 
            end if 
        end do 


    end function

    !------------------------------------------------------------------!
    !                          Auxiliary functions                     !
    !------------------------------------------------------------------!

    ! Box checker 
    function IsInBox(minx, maxx, miny, maxy, x, y) result(in) 

        ! Description
        !============
        ! Simple box checker. x and by values of the box should be given
        ! in maxx, maxy, ... , values should be dimension (:).
        ! Results may vary when values are up to precision on the box, 
        ! since this precision is not accounted for.

        ! Declare variables
        !==================
        ! Arguments
        real(R8)                    :: minx, maxx, miny, maxy
        real(R8), allocatable, intent(in)       :: x(:), y(:)
        logical, allocatable        :: in(:)

        ! Auxiliary

        ! Initialize
        !===========
        ! Allocate
        allocate(in(size(x)))

        ! Check if in box
        !================
        in = ( (x >= minx) .and. (x <= maxx ) ) .and. &
            ( (y >= miny) .and. (y <= maxy) )

    end function 

end module
