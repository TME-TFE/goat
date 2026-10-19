!======================================================================!
!                                                                      !
!                            DOCUMENTATION                             !
!                                                                      !
!======================================================================!

! Description
!============
! This module reworks the volume regions (cvReg) and face regions
! (fcReg) of a Goat grid into the region numbering that
! SOLPS-ITER expects, based on the catalog topology recognized by
! ClassifySOLPSCatalogTopology.
! One routine per catalog family. Currently implemented:
!  - SetSOLPSRegionsLIM: limiter topologies (NREG = 2, FREG = 6)
!  - SetSOLPSRegionsSN: lower/upper single null (NREG = 4, FREG = 13)
!  - SetSOLPSRegionsDN: CDN and DDN bottom/top (NREG = 8,
!    FREG = 26/27)

module ggmod_solpsregions

    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_definitions
    use mod_constants, only: posinfval_R8
    use goatmod_types, only: GridUDT, VesselUDT, GetVertFace, &
        GetFaceCell, GetCellFace
    use ggmod_topology2D

    implicit none
    private
    public :: SetSOLPSRegionsLIM, SetSOLPSRegionsSN, &
        SetSOLPSRegionsDN, AddSecondaryXPointRegions

contains

    ! SOLPS region conventions for the LIM (limiter) catalog topology
    subroutine SetSOLPSRegionsLIM(simgrid, topomesh, vessel)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! limiter case (NREG = 2, FREG = 6):
        !
        ! Volume regions (cvReg):
        !   1 = core (inside the separatrix)
        !   2 = SOL (everything else)
        !
        ! Face regions (fcReg):
        !   1 = left target: the part of the touched wall structure
        !       on the clockwise side (standard (R,z) view) of the
        !       point where the separatrix touches the wall.
        !   2 = right target: the counterclockwise part.
        !   3 = (empty) core connection.
        !   4 = core boundary (the inner grid boundary)
        !   5 = separatrix: internal aligned faces on the interface
        !       between core and SOL topomesh cells
        !   6 = the rest of the wall / outer boundary
        !
        ! The targets are bounded to the composite limiter surface
        ! formed by the SOLPS limiter target segments
        ! (goat.vessel.solps_limiter_target): the left target extends
        ! from the separatrix-wall touch point to one end of that
        ! composite surface, the right target to the other end. For a
        ! SOLPS limiter grid that list is REQUIRED to be non-empty
        ! and to form exactly one consecutive stretch of wall
        ! elements - otherwise this routine errors out (the same
        ! requirement is already enforced fail-fast, before any
        ! gridding, in InitializeGGTMLineRefiner).

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel

        ! Auxiliary
        integer(I8), allocatable        :: cc(:), ccf(:), tp(:), &
            tcv(:), tfc(:), tf2(:)
        logical, allocatable            :: coremaskTM(:), &
            iscorefTM(:), issepTM(:), slmask(:)
        integer(I8)                     :: i, k, itp, istrtouch, lbl, &
            nleft, nright, nwall, ncoreb, nsep, nmix, iv, nrun, kp
        real(R8)                        :: xc, yc, xt, yt, xm, ym, &
            d, dmin, tht, thm, dth
        logical                         :: havetargets

        ! Volume regions
        !===============
        ! Core (marked with the SOLPS core region IDs during the
        ! label translation, or already 1/2 from the generic limiter
        ! remapping - the test is idempotent) -> 1, rest -> 2
        where (mod(simgrid%cell%reg - SOLPScoreregID, &
            SOLPScoreregIDincr) == 0)
            simgrid%cell%reg = 1
        elsewhere
            simgrid%cell%reg = 2
        end where

        ! Topomesh masks
        !===============
        allocate(coremaskTM(topomesh%cell%ntot), &
            iscorefTM(topomesh%face%ntot), &
            issepTM(topomesh%face%ntot))
        coremaskTM = .false.
        iscorefTM = .false.
        issepTM = .false.
        cc = topomesh%GetCoreCellIDs()
        coremaskTM(cc) = .true.
        ccf = topomesh%GetCoreFaceIDs()
        if (size(ccf) > 0) iscorefTM(ccf) = .true.

        ! Separatrix topomesh faces: interfaces between a core and a
        ! non-core topomesh cell
        do i = 1, topomesh%face%ntot
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (coremaskTM(tfc(1)) .neqv. coremaskTM(tfc(2))) then
                    issepTM(i) = .true.
                end if
            end if
        end do

        ! Core centroid (reference for the poloidal angle)
        !==================================================
        xc = 0.0_R8
        yc = 0.0_R8
        k = 0
        do i = 1, size(cc)
            tcv = topomesh%cell%GetVert(cc(i))
            xc = xc + sum(topomesh%vert%x(tcv))
            yc = yc + sum(topomesh%vert%y(tcv))
            k = k + size(tcv)
        end do
        xc = xc/real(k, kind=R8)
        yc = yc/real(k, kind=R8)

        ! Separatrix-wall touch point and touched wall structure
        !========================================================
        ! The touch point is the topomesh vertex shared by a
        ! separatrix (core/SOL interface) face and a wall (boundary)
        ! face. This is robust against the tangency-point typing and
        ! closed-contour checks of the topomesh construction.
        havetargets = .false.
        xt = 0.0_R8
        yt = 0.0_R8
        istrtouch = 0

        ! Requirement for SOLPS limiter grids: the SOLPS limiter
        ! target segments (goat.vessel.solps_limiter_target) must be
        ! given and must form exactly one consecutive stretch of wall
        ! elements (checked modulo the wall wrap-around)
        allocate(slmask(vessel%nstructures))
        slmask = .false.
        if (vessel%nsolpslim <= 0) then
            call gdErrorHandler('SetSOLPSRegionsLIM: ' // &
                'goat.vessel.solps_limiter_target is empty - a ' // &
                'SOLPS limiter grid requires the limiter ' // &
                'surface to be given there as one consecutive ' // &
                'stretch of wall elements')
        end if
        do i = 1, int(vessel%nsolpslim, kind=I8)
            k = int(vessel%solpslimind(i), kind=I8)
            if (slmask(k)) then
                call gdErrorHandler('SetSOLPSRegionsLIM: ' // &
                    'duplicate wall element in ' // &
                    'goat.vessel.solps_limiter_target')
            end if
            slmask(k) = .true.
        end do
        nrun = 0
        do k = 1, int(vessel%nstructures, kind=I8)
            kp = k - 1
            if (kp == 0) kp = int(vessel%nstructures, kind=I8)
            if (slmask(k) .and. .not. slmask(kp)) nrun = nrun + 1
        end do
        if ((nrun /= 1) .and. &
            (count(slmask) /= int(vessel%nstructures, kind=I8))) &
            then
            call gdErrorHandler('SetSOLPSRegionsLIM: the wall ' // &
                'elements in goat.vessel.solps_limiter_target ' // &
                'must form exactly one consecutive stretch ' // &
                '(the limiter surface) for a SOLPS limiter grid')
        end if

        allocate(tp(0))
        do i = 1, topomesh%face%ntot
            if (.not. issepTM(i)) cycle
            do k = 1, 2
                iv = topomesh%face%vert(i, k)
                if (size(tp) > 0) then
                    if (any(tp == iv)) cycle
                end if
                tf2 = topomesh%vert%GetFace(iv)
                if (any((topomesh%face%type(tf2) == TMfacebndID) &
                        .or. &
                        (topomesh%face%type(tf2) == TMfacealbndID))) &
                    then
                    tp = [tp, iv]
                end if
            end do
        end do
        if (size(tp) == 0) then
            print *, 'SetSOLPSRegionsLIM: no separatrix-wall ' // &
                'touch point found - leaving targets empty'
        else
            if (size(tp) > 1) then
                print *, 'SetSOLPSRegionsLIM: ', size(tp), &
                    ' separatrix-wall touch vertices found, ' // &
                    'using the first: ', tp
            end if
            itp = tp(1)
            xt = topomesh%vert%x(itp)
            yt = topomesh%vert%y(itp)

            ! Touched structure: nearest structure.dat sub-structure
            dmin = posinfval_R8()
            do k = 1, int(vessel%nstructures, kind=I8)
                d = DistToStructure(xt, yt, k)
                if (d < dmin) then
                    dmin = d
                    istrtouch = k
                end if
            end do
            havetargets = .true.
            print *, 'SetSOLPSRegionsLIM: separatrix touch point ', &
                xt, yt, ' on structure ', istrtouch
            if (.not. slmask(istrtouch)) then
                print *, 'SetSOLPSRegionsLIM: WARNING: the touch ' // &
                    'point does not lie on the limiter surface ' // &
                    '(goat.vessel.solps_limiter_target) - one of ' // &
                    'the two targets will be empty'
            end if
        end if

        ! Face regions
        !=============
        tht = atan2(yt - yc, xt - xc)
        simgrid%face%reg = 0
        nleft = 0
        nright = 0
        nwall = 0
        ncoreb = 0
        nsep = 0
        do i = 1, simgrid%face%ntot
            lbl = simgrid%face%TMfacelabel(i)
            if (simgrid%face%BF(i)) then

                ! Core boundary
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (iscorefTM(lbl)) then
                        simgrid%face%reg(i) = 4
                        ncoreb = ncoreb + 1
                        cycle
                    end if
                end if

                ! Wall: split the touched structure into left/right
                ! targets, everything else is region 6
                simgrid%face%reg(i) = 6
                nwall = nwall + 1
                if (havetargets) then
                    xm = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%x(simgrid%face%vert(i, 2)))
                    ym = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%y(simgrid%face%vert(i, 2)))
                    dmin = posinfval_R8()
                    k = 0
                    do lbl = 1, int(vessel%nstructures, kind=I8)
                        d = DistToStructure(xm, ym, lbl)
                        if (d < dmin) then
                            dmin = d
                            k = lbl
                        end if
                    end do
                    if (slmask(k)) then
                        ! Face on the composite limiter surface:
                        ! clockwise side of the touch point = left
                        ! target (region 1), counterclockwise =
                        ! right target (region 2)
                        thm = atan2(ym - yc, xm - xc)
                        dth = atan2(sin(thm - tht), cos(thm - tht))
                        nwall = nwall - 1
                        if (dth < 0.0_R8) then
                            simgrid%face%reg(i) = 1
                            nleft = nleft + 1
                        else
                            simgrid%face%reg(i) = 2
                            nright = nright + 1
                        end if
                    end if
                end if

            else

                ! Internal faces: separatrix -> 5, rest stays 0
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (issepTM(lbl)) then
                        simgrid%face%reg(i) = 5
                        nsep = nsep + 1
                    end if
                end if

            end if
        end do

        print *, 'SetSOLPSRegionsLIM: fcReg counts: left target ', &
            nleft, ', right target ', nright, ', core boundary ', &
            ncoreb, ', separatrix ', nsep, ', wall/outer ', nwall

        ! Flux tube regions
        !==================
        ! In a limiter topology every flux tube lies entirely in the
        ! core or entirely in the SOL: inherit the (rewritten) cell
        ! region -> 1 = core, 2 = SOL. No other values exist.
        nmix = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            k = fd%fluxtubecells(fd%fluxtubecellsP(i, 1))
            fd%fluxtuberegID(i) = simgrid%cell%reg(k)
            do lbl = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (simgrid%cell%reg(fd%fluxtubecells(lbl)) /= &
                    fd%fluxtuberegID(i)) nmix = nmix + 1
            end do
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsLIM: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL regions ' // &
                '(tube crossing the separatrix?)'
        end if
        print *, 'SetSOLPSRegionsLIM: ftReg: ', &
            count(fd%fluxtuberegID(1:fd%nFt) == 1), ' core tubes, ', &
            count(fd%fluxtuberegID(1:fd%nFt) == 2), ' SOL tubes'
        end associate

    contains

        ! Minimum distance from a point to vessel sub-structure ks
        function DistToStructure(px, py, ks) result(dres)
            real(R8), intent(in)        :: px, py
            integer(I8), intent(in)     :: ks
            real(R8)                    :: dres
            integer(I8)                 :: jj, npt
            real(R8)                    :: x1, y1, x2, y2
            dres = posinfval_R8()
            npt = vessel%structures(ks)%np
            do jj = 1, npt - 1
                x1 = vessel%structures(ks)%x(jj)
                y1 = vessel%structures(ks)%y(jj)
                x2 = vessel%structures(ks)%x(jj + 1)
                y2 = vessel%structures(ks)%y(jj + 1)
                dres = min(dres, DistToSegment(px, py, x1, y1, x2, y2))
            end do
            if (vessel%structures(ks)%isclosed .and. (npt > 2)) then
                x1 = vessel%structures(ks)%x(npt)
                y1 = vessel%structures(ks)%y(npt)
                x2 = vessel%structures(ks)%x(1)
                y2 = vessel%structures(ks)%y(1)
                dres = min(dres, DistToSegment(px, py, x1, y1, x2, y2))
            end if
        end function

        ! Minimum distance from a point to one segment
        function DistToSegment(px, py, x1, y1, x2, y2) result(dres)
            real(R8), intent(in)    :: px, py, x1, y1, x2, y2
            real(R8)                :: dres
            real(R8)                :: dx, dy, t
            dx = x2 - x1
            dy = y2 - y1
            if ((dx*dx + dy*dy) <= 0.0_R8) then
                dres = sqrt((px - x1)**2 + (py - y1)**2)
                return
            end if
            t = ((px - x1)*dx + (py - y1)*dy)/(dx*dx + dy*dy)
            t = max(0.0_R8, min(1.0_R8, t))
            dres = sqrt((px - (x1 + t*dx))**2 + (py - (y1 + t*dy))**2)
        end function

    end subroutine

    ! SOLPS region conventions for the SN (single null) catalog
    ! topologies
    subroutine SetSOLPSRegionsSN(simgrid, topomesh, vessel, isupper)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! single-null case (NREG = 4, FREG = 13). West/East follow the
        ! convention with its reversed target ordering: West = inboard
        ! for a lower single null, outboard for an upper single null.
        !
        ! The construction starts from the 8 "rays" leaving the
        ! X-point grid vertex: 4 are the separatrix branches (aligned
        ! faces), the other 4 are chains of non-aligned (radial) grid
        ! faces walked outwards from the X-point until a boundary
        ! vertex: into the core (core cut, face region 5), into the
        ! PFR (divertor cut, face region 6, splitting the inner and
        ! outer divertor), and into the two SOL sectors (divertor
        ! entrances, face regions 2/3, splitting main SOL from
        ! divertor SOL).
        !
        ! Volume regions (cvReg):
        !   1 = core
        !   2 = main SOL
        !   3 = Western divertor (divertor SOL + PFR, West side)
        !   4 = Eastern divertor (divertor SOL + PFR, East side)
        !
        ! Face regions (fcReg):
        !   1 = Western target: wall faces on the structures given in
        !       goat.vessel.solps_lower_divertor_inner_target (lower
        !       SN) / solps_upper_divertor_outer_target (upper SN).
        !       The West and East target lists may share a structure
        !       (both legs striking one flat wall run): its faces are
        !       then split between regions 1 and 4 at the divertor
        !       cut (via the volume regions), and the PFR wall
        !       regions 7/9 in between do not exist.
        !   2 = entrance to the Western divertor (radial chain from
        !       the X-point to the wall through the West SOL)
        !   3 = entrance to the Eastern divertor
        !   4 = Eastern target (solps_lower_divertor_outer_target /
        !       solps_upper_divertor_inner_target)
        !   5 = core cut (radial chain X-point -> core boundary)
        !   6 = divertor cut (radial chain X-point -> wall through
        !       the PFR, splitting the two divertors)
        !   7 = Western divertor PFR wall: boundary faces adjacent to
        !       PFR cells of region 3 that are not target faces
        !   8 = core boundary (boundary faces adjacent to core cells)
        !   9 = Eastern divertor PFR wall (as 7, for region 4)
        !  10 = separatrix (internal faces on the core-enclosing
        !       part of the separatrix only - the divertor legs
        !       are excluded)
        !  11 = Western divertor main-chamber wall: boundary faces
        !       adjacent to divertor-SOL cells of region 3 that are
        !       not target faces (empty when the entrance chain
        !       lands on the target itself)
        !  12 = main-chamber (main SOL) wall and any remaining wall
        !  13 = Eastern divertor main-chamber wall (as 11, region 4)
        !
        ! Flux tube regions (ftReg): 1 = core, 2 = SOL, 3 = PFR.
        !
        ! Empty divertor target inputs are tolerated (warning): the
        ! corresponding target face region stays empty and its wall
        ! faces fall to regions 7/9/11/13/12. Faces are attributed to
        ! a wall structure only when they lie within dtolwall of it,
        ! so flux-aligned domain boundaries far from the wall cannot
        ! be mistaken for targets.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel
        logical, intent(in)             :: isupper

        ! Parameters: cell kinds and wall attribution tolerance
        integer(I8), parameter          :: KCORE = 1, KSOL = 2, &
            KPFR = 3
        real(R8), parameter             :: dtolwall = 5.0e-3_R8

        ! Auxiliary
        integer(I8), allocatable        :: cc(:), xp(:), tcv(:), &
            tfc(:), vf(:), cf(:), cvkind(:), comp(:), queue(:), &
            chainreg(:), compreg(:), compn(:), zn(:), zkind(:)
        logical, allocatable            :: issepTM(:), ismain(:), &
            istarget1(:), istarget4(:), ispfrhalf(:)
        real(R8), allocatable           :: compx(:), compy(:), &
            zpsi(:), zcx(:), zcy(:)
        integer(I8)                     :: i, j, k, ich, iv, ivX, w, &
            f, fprev, lbl, c, c2, ncomp, nsolcomp, nqh, nqt, nmix, &
            nsepal, npfrhalf, fstart(4), chlen(4), termvx(4), &
            ncand(4), nfreg(0:13), ncvreg(4), nftreg(3), nprob1, iprim
        real(R8)                        :: xc, yc, ux, uy, xxp, yxp, &
            d, dmin, score, dirx, diry, dp, s, psiX, psicore, &
            xm, ym, vx, vy, cxr, cyr, ycref, fcref
        logical                         :: okref
        character(8)                    :: chname(4)

        chname = ['core    ', 'PFR     ', 'West SOL', 'East SOL']
        print *, 'SetSOLPSRegionsSN: applying SOLPS single-null ' // &
            'regions (orientation: ' // &
            trim(merge('upper', 'lower', isupper)) // ')'

        ! X-point and reference directions
        !=================================
        xp = topomesh%GetXPointIDs()
        ! Narrow-grid rule: an X-point whose flux surface was removed
        ! from the delivered grid by the narrow-grid target filter
        ! (orphaned grid vertex, fieldlineID 0) does not exist for the
        ! SOLPS region conventions - drop it from the working list. On
        ! a complete (wide) grid every X-point stays.
        if (size(xp) > 0) then
            if (any(simgrid%vert%fieldlineID(xp) <= 0)) then
                print *, 'SetSOLPSRegionsSN: ignoring ', &
                    count(simgrid%vert%fieldlineID(xp) <= 0), &
                    ' X-point(s) without a flux surface in the ' // &
                    'filtered grid'
                xp = pack(xp, simgrid%vert%fieldlineID(xp) > 0)
            end if
        end if
        if (size(xp) == 0) then
            print *, 'SetSOLPSRegionsSN: no X-point in the ' // &
                'topomesh - leaving the generic regions in place'
            return
        end if
        ! Primary X-point = the innermost on the single-null side; any
        ! further same-side X-points (snowflake companions) are handled
        ! afterwards by AddSecondaryXPointRegions.
        call SOLPSReferenceCenterFlux(topomesh, ycref, fcref, okref)
        if (.not. okref) then
            print *, 'SetSOLPSRegionsSN: no core reference - ' // &
                'leaving the generic regions in place'
            return
        end if
        iprim = 0
        do i = 1, size(xp)
            if ((topomesh%vert%y(xp(i)) >= ycref) .neqv. isupper) cycle
            d = abs(topomesh%vert%fval(xp(i)) - fcref)
            if ((iprim == 0) .or. (d < dmin)) then
                iprim = i
                dmin = d
            end if
        end do
        if (iprim == 0) then
            print *, 'SetSOLPSRegionsSN: no X-point on the ' // &
                trim(merge('upper', 'lower', isupper)) // ' side - ' // &
                'leaving the generic regions in place'
            return
        end if
        if (size(xp) > 1) then
            print *, 'SetSOLPSRegionsSN: ', size(xp), ' X-points; ' // &
                'primary = innermost on the ' // &
                trim(merge('upper', 'lower', isupper)) // ' side ' // &
                '(companions handled by AddSecondaryXPointRegions)'
        end if
        xxp = topomesh%vert%x(xp(iprim))
        yxp = topomesh%vert%y(xp(iprim))

        ! Grid vertex sitting on the X-point
        ivX = 1
        dmin = posinfval_R8()
        do i = 1, simgrid%vert%ntot
            d = (simgrid%vert%x(i) - xxp)**2 + &
                (simgrid%vert%y(i) - yxp)**2
            if (d < dmin) then
                dmin = d
                ivX = i
            end if
        end do
        if (sqrt(dmin) > 1.0e-6_R8) then
            print *, 'SetSOLPSRegionsSN: WARNING: nearest grid ' // &
                'vertex to the X-point is ', sqrt(dmin), ' m away'
        end if
        psiX = simgrid%vert%psi(ivX)

        ! Fail-fast: the primary X-point's strike legs must land on the
        ! declared divertor targets (both target lists must be given).
        ! West/East follow the SN convention (West = inner for a lower
        ! SN, outer for an upper SN), but the guard only needs the two
        ! lists of the X-point's side.
        ! Check every X-point on this single-null side (primary + any
        ! snowflake companions); the guard keeps only the real divertor
        ! legs (two per X-point).
        if (isupper) then
            call GuardPrimaryStrikes(topomesh, vessel, dtolwall, &
                'SetSOLPSRegionsSN', xp, [3_I8, 4_I8], nprob1)
        else
            call GuardPrimaryStrikes(topomesh, vessel, dtolwall, &
                'SetSOLPSRegionsSN', xp, [1_I8, 2_I8], nprob1)
        end if
        if (nprob1 > 0) then
            call gdErrorHandler('SetSOLPSRegionsSN: primary strike ' // &
                'points are not all on declared divertor targets ' // &
                '(see the messages above) - set/split the ' // &
                'goat.vessel.solps_{lower,upper}_divertor_' // &
                '{inner,outer}_target lists so every primary strike ' // &
                'lands on target wall')
        end if

        ! Core centroid: defines the reference direction
        ! u = (X-point -> core), whose left side is SOLPS West
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            print *, 'SetSOLPSRegionsSN: no core topomesh cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        xc = 0.0_R8
        yc = 0.0_R8
        k = 0
        do i = 1, size(cc)
            tcv = topomesh%cell%GetVert(cc(i))
            xc = xc + sum(topomesh%vert%x(tcv))
            yc = yc + sum(topomesh%vert%y(tcv))
            k = k + size(tcv)
        end do
        xc = xc/real(k, kind=R8)
        yc = yc/real(k, kind=R8)
        ux = xc - xxp
        uy = yc - yxp

        ! Separatrix topomesh faces: face region 10
        ! is only the core-enclosing part of the
        ! separatrix, so the interfaces between a core and
        ! a non-core topomesh cell - the divertor legs
        ! (PFR/SOL interfaces) are excluded
        allocate(issepTM(topomesh%face%ntot))
        issepTM = .false.
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= TMfacesepID) cycle
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (any(tfc(1) == cc) .neqv. any(tfc(2) == cc)) then
                    issepTM(i) = .true.
                end if
            end if
        end do

        ! Cell kinds
        !===========
        ! Core from the region marks set during the label
        ! translation; PFR/SOL per ZONE (the pre-overwrite cvReg
        ! values, i.e. the topomesh cells, which never straddle the
        ! separatrix) from the mean flux side relative to the
        ! X-point: the PFR flux is always core-side of its null.
        ! The zone average makes the test immune to interpolation
        ! noise on the cells hugging the separatrix.
        allocate(zpsi(maxval(simgrid%cell%reg)), &
            zn(maxval(simgrid%cell%reg)), &
            zkind(maxval(simgrid%cell%reg)), &
            zcx(maxval(simgrid%cell%reg)), &
            zcy(maxval(simgrid%cell%reg)))
        zpsi = 0.0_R8
        zcx = 0.0_R8
        zcy = 0.0_R8
        zn = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            zpsi(k) = zpsi(k) + simgrid%cell%psi(i)
            zcx(k) = zcx(k) + simgrid%cell%x(i)
            zcy(k) = zcy(k) + simgrid%cell%y(i)
            zn(k) = zn(k) + 1
        end do
        psicore = 0.0_R8
        k = 0
        do j = 1, size(zkind)
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                psicore = psicore + zpsi(j)
                k = k + zn(j)
            end if
        end do
        if (k == 0) then
            print *, 'SetSOLPSRegionsSN: no core grid cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        s = sign(1.0_R8, psiX - psicore/real(k, kind=R8))
        ! PFR = flux core-side of the null AND on the divertor side of
        ! the X-point (behind it, away from the core). The second
        ! condition is redundant for a genuine single null, but on a
        ! reclassified snowflake (narrow grid whose second divertor was
        ! filtered away) the surviving companion separatrix runs a
        ! fraction of a permille outside the primary one: upstream
        ! zones of that inter-separatrix band have mean flux within
        ! interpolation noise of psiX, and a flux-only test scatters
        ! them into the PFR kind - they then flood-connect through the
        ! X-region into a divertor region (cells of the main SOL ring
        ! labeled as divertor all around the machine)
        do j = 1, size(zkind)
            zkind(j) = KSOL
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                zkind(j) = KCORE
            elseif ((s*(zpsi(j)/real(zn(j), kind=R8) - psiX) < &
                0.0_R8) .and. &
                (((zcx(j)/real(zn(j), kind=R8) - xxp)*ux + &
                  (zcy(j)/real(zn(j), kind=R8) - yxp)*uy) < &
                0.0_R8)) then
                zkind(j) = KPFR
            end if
        end do
        allocate(cvkind(simgrid%cell%ntot))
        do i = 1, simgrid%cell%ntot
            cvkind(i) = zkind(simgrid%cell%reg(i))
        end do

        ! The four radial rays leaving the X-point vertex
        !================================================
        ! Non-aligned internal faces at the X-point vertex, classified
        ! by the sector they walk into: core/PFR by the flux side
        ! (told apart by the dot product with u), West/East SOL by
        ! the cross product with u
        fstart = 0
        ncand = 0
        nsepal = 0
        vf = GetVertFace(simgrid%vert, ivX)
        do j = 1, size(vf)
            f = vf(j)
            if (simgrid%face%BF(f)) cycle
            if (simgrid%face%aligned(f) == 1) then
                nsepal = nsepal + 1
                cycle
            end if
            if (simgrid%face%vert(f, 1) == ivX) then
                w = simgrid%face%vert(f, 2)
            else
                w = simgrid%face%vert(f, 1)
            end if
            dirx = simgrid%vert%x(w) - xxp
            diry = simgrid%vert%y(w) - yxp
            dp = s*(simgrid%vert%psi(w) - psiX)
            if (dp < 0.0_R8) then
                if ((dirx*ux + diry*uy) > 0.0_R8) then
                    ich = 1
                else
                    ich = 2
                end if
            else
                if ((ux*diry - uy*dirx) > 0.0_R8) then
                    ich = 3
                else
                    ich = 4
                end if
            end if
            ncand(ich) = ncand(ich) + 1
            if (fstart(ich) == 0) fstart(ich) = f
        end do
        print *, 'SetSOLPSRegionsSN: X-point grid vertex ', ivX, &
            ' at ', xxp, yxp, ': ', nsepal, ' separatrix rays, ' // &
            'radial ray candidates (core/PFR/West/East): ', ncand
        do ich = 1, 4
            if (fstart(ich) == 0) then
                print *, 'SetSOLPSRegionsSN: WARNING: no radial ' // &
                    'ray into the ' // trim(chname(ich)) // &
                    ' sector at the X-point'
            end if
        end do

        ! Walk the four radial chains from the X-point to a boundary
        ! vertex. Their faces are tagged with the face region of the
        ! cut they realize: 5 = core cut, 6 = divertor (PFR) cut,
        ! 2 = West divertor entrance, 3 = East divertor entrance
        allocate(chainreg(simgrid%face%ntot))
        chainreg = 0
        chlen = 0
        termvx = 0
        do ich = 1, 4
            if (fstart(ich) == 0) cycle
            iv = ivX
            f = fstart(ich)
            do
                select case (ich)
                case (1)
                    chainreg(f) = 5
                case (2)
                    chainreg(f) = 6
                case (3)
                    chainreg(f) = 2
                case (4)
                    chainreg(f) = 3
                end select
                chlen(ich) = chlen(ich) + 1
                if (chlen(ich) > simgrid%face%ntot) then
                    print *, 'SetSOLPSRegionsSN: WARNING: runaway ' // &
                        trim(chname(ich)) // ' chain walk - aborted'
                    exit
                end if
                if (simgrid%face%vert(f, 1) == iv) then
                    w = simgrid%face%vert(f, 2)
                else
                    w = simgrid%face%vert(f, 1)
                end if
                if (simgrid%vert%BV(w)) then
                    termvx(ich) = w
                    exit
                end if

                ! Continue on the unvisited radial face that best
                ! keeps the walking direction
                vx = simgrid%vert%x(w) - simgrid%vert%x(iv)
                vy = simgrid%vert%y(w) - simgrid%vert%y(iv)
                fprev = f
                f = 0
                score = -posinfval_R8()
                cf = GetVertFace(simgrid%vert, w)
                do j = 1, size(cf)
                    if (cf(j) == fprev) cycle
                    if (simgrid%face%BF(cf(j))) cycle
                    if (simgrid%face%aligned(cf(j)) == 1) cycle
                    if (chainreg(cf(j)) /= 0) cycle
                    if (simgrid%face%vert(cf(j), 1) == w) then
                        c2 = simgrid%face%vert(cf(j), 2)
                    else
                        c2 = simgrid%face%vert(cf(j), 1)
                    end if
                    dirx = simgrid%vert%x(c2) - simgrid%vert%x(w)
                    diry = simgrid%vert%y(c2) - simgrid%vert%y(w)
                    d = (dirx*vx + diry*vy)/ &
                        max(sqrt(dirx**2 + diry**2), tiny(1.0_R8))
                    if (d > score) then
                        score = d
                        f = cf(j)
                    end if
                end do
                if (f == 0) then
                    print *, 'SetSOLPSRegionsSN: WARNING: ' // &
                        trim(chname(ich)) // ' chain stopped at ' // &
                        'internal vertex ', w, ' after ', &
                        chlen(ich), ' faces'
                    termvx(ich) = w
                    exit
                end if
                iv = w
            end do
        end do
        do ich = 1, 4
            if (termvx(ich) == 0) cycle
            xm = simgrid%vert%x(termvx(ich))
            ym = simgrid%vert%y(termvx(ich))
            k = 0
            dmin = posinfval_R8()
            do j = 1, int(vessel%nstructures, kind=I8)
                d = StructureDistance(vessel, xm, ym, j)
                if (d < dmin) then
                    dmin = d
                    k = j
                end if
            end do
            if (ich == 1) then
                print *, 'SetSOLPSRegionsSN: ', trim(chname(ich)), &
                    ' cut: ', chlen(ich), ' faces, ends at ', xm, ym
            else
                print *, 'SetSOLPSRegionsSN: ', trim(chname(ich)), &
                    ' cut: ', chlen(ich), ' faces, intercepts ' // &
                    'the wall at ', xm, ym, ' (structure ', k, ')'
            end if
        end do

        ! Volume regions
        !===============
        ! Flood fill of the SOL cells with the two entrance chains as
        ! barriers (main SOL + the two divertor SOLs), then of the
        ! PFR cells with the divertor cut as barrier (the two PFR
        ! halves)
        allocate(comp(simgrid%cell%ntot), queue(simgrid%cell%ntot))
        comp = 0
        ncomp = 0
        call FloodFill(KSOL, 2_I8, 3_I8)
        nsolcomp = ncomp
        call FloodFill(KPFR, 6_I8, 6_I8)
        if (nsolcomp /= 3) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nsolcomp, &
                ' SOL components instead of 3 (incomplete ' // &
                'entrance chains?)'
        end if

        allocate(ismain(ncomp), compreg(ncomp), compx(ncomp), &
            compy(ncomp), compn(ncomp), ispfrhalf(ncomp))
        ismain = .false.
        compreg = 0

        ! The two PFR halves are the components adjacent to the
        ! divertor cut. Any other flux-side component is a far-field
        ! pocket whose flux folded back to the core side (e.g.
        ! beyond a clipped secondary separatrix): those belong to
        ! the main SOL, not to a divertor
        ispfrhalf = .false.
        do i = 1, simgrid%face%ntot
            if (chainreg(i) /= 6) cycle
            tfc = GetFaceCell(simgrid%face, i)
            do k = 1, size(tfc)
                if (cvkind(tfc(k)) == KPFR) then
                    ispfrhalf(comp(tfc(k))) = .true.
                end if
            end do
        end do
        npfrhalf = 0
        do c = nsolcomp + 1, ncomp
            if (ispfrhalf(c)) npfrhalf = npfrhalf + 1
        end do
        nmix = 0
        do i = 1, simgrid%cell%ntot
            if (cvkind(i) /= KPFR) cycle
            if (ispfrhalf(comp(i))) cycle
            cvkind(i) = KSOL
            ismain(comp(i)) = .true.
            nmix = nmix + 1
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: ', nmix, ' flux-side ' // &
                'cells in ', ncomp - nsolcomp - npfrhalf, &
                ' components not connected to the PFR (folded ' // &
                'far-field flux?) - assigned to the main SOL'
        end if
        if (npfrhalf /= 2) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', npfrhalf, &
                ' PFR halves instead of 2 (incomplete divertor ' // &
                'cut?) - splitting the PFR by the side of the ' // &
                'X-point -> core axis instead'
        end if

        ! Main SOL = the SOL component(s) radially adjacent to the
        ! core across the separatrix; the divertor SOL components
        ! are identified by the entrance chain bounding them
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) cycle
            tfc = GetFaceCell(simgrid%face, i)
            if (size(tfc) /= 2) cycle
            if ((cvkind(tfc(1)) == KCORE) .and. &
                (cvkind(tfc(2)) == KSOL)) then
                ismain(comp(tfc(2))) = .true.
            elseif ((cvkind(tfc(2)) == KCORE) .and. &
                (cvkind(tfc(1)) == KSOL)) then
                ismain(comp(tfc(1))) = .true.
            end if
        end do
        do i = 1, simgrid%face%ntot
            if ((chainreg(i) /= 2) .and. (chainreg(i) /= 3)) cycle
            tfc = GetFaceCell(simgrid%face, i)
            do k = 1, size(tfc)
                if (cvkind(tfc(k)) /= KSOL) cycle
                c = comp(tfc(k))
                if (ismain(c)) cycle
                if (chainreg(i) == 2) then
                    compreg(c) = 3
                else
                    compreg(c) = 4
                end if
            end do
        end do

        ! PFR components: West (region 3) / East (region 4) by the
        ! side of the X-point -> core direction
        compx = 0.0_R8
        compy = 0.0_R8
        compn = 0
        do i = 1, simgrid%cell%ntot
            if (cvkind(i) /= KPFR) cycle
            compx(comp(i)) = compx(comp(i)) + simgrid%cell%x(i)
            compy(comp(i)) = compy(comp(i)) + simgrid%cell%y(i)
            compn(comp(i)) = compn(comp(i)) + 1
        end do
        do c = nsolcomp + 1, ncomp
            if (.not. ispfrhalf(c)) cycle
            if (compn(c) == 0) cycle
            cxr = compx(c)/real(compn(c), kind=R8) - xxp
            cyr = compy(c)/real(compn(c), kind=R8) - yxp
            if ((ux*cyr - uy*cxr) > 0.0_R8) then
                compreg(c) = 3
            else
                compreg(c) = 4
            end if
        end do

        ! Cell region assignment
        nmix = 0
        do i = 1, simgrid%cell%ntot
            select case (cvkind(i))
            case (KCORE)
                simgrid%cell%reg(i) = 1
            case (KSOL)
                if (ismain(comp(i))) then
                    simgrid%cell%reg(i) = 2
                elseif (compreg(comp(i)) /= 0) then
                    simgrid%cell%reg(i) = compreg(comp(i))
                else
                    simgrid%cell%reg(i) = 2
                    nmix = nmix + 1
                end if
            case default
                if ((npfrhalf == 2) .and. &
                    (compreg(comp(i)) /= 0)) then
                    simgrid%cell%reg(i) = compreg(comp(i))
                else
                    ! Fallback split by the X-point -> core axis
                    cxr = simgrid%cell%x(i) - xxp
                    cyr = simgrid%cell%y(i) - yxp
                    if ((ux*cyr - uy*cxr) > 0.0_R8) then
                        simgrid%cell%reg(i) = 3
                    else
                        simgrid%cell%reg(i) = 4
                    end if
                end if
            end select
        end do
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nmix, &
                ' SOL cells in components touching no entrance ' // &
                'chain and no core - assigned to the main SOL'
        end if
        ncvreg = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k >= 1) .and. (k <= 4)) ncvreg(k) = ncvreg(k) + 1
        end do
        print *, 'SetSOLPSRegionsSN: cvReg counts: core ', &
            ncvreg(1), ', main SOL ', ncvreg(2), ', West ' // &
            'divertor ', ncvreg(3), ', East divertor ', ncvreg(4)

        ! Face regions
        !=============
        ! Western/Eastern target structures from the Goat inputs,
        ! with reversed ordering: lower SN -> West =
        ! inner target, upper SN -> West = outer target
        allocate(istarget1(vessel%nstructures), &
            istarget4(vessel%nstructures))
        istarget1 = .false.
        istarget4 = .false.
        if (isupper) then
            do i = 1, int(vessel%nsolpsudo, kind=I8)
                istarget1(vessel%solpsudoind(i)) = .true.
            end do
            do i = 1, int(vessel%nsolpsudi, kind=I8)
                istarget4(vessel%solpsudiind(i)) = .true.
            end do
        else
            do i = 1, int(vessel%nsolpsldi, kind=I8)
                istarget1(vessel%solpsldiind(i)) = .true.
            end do
            do i = 1, int(vessel%nsolpsldo, kind=I8)
                istarget4(vessel%solpsldoind(i)) = .true.
            end do
        end if
        if (count(istarget1) == 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: no Western ' // &
                'target structures given (goat.vessel.solps_' // &
                trim(merge('upper_divertor_outer', &
                'lower_divertor_inner', isupper)) // '_target) - ' // &
                'face region 1 will be empty'
        end if
        if (count(istarget4) == 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: no Eastern ' // &
                'target structures given (goat.vessel.solps_' // &
                trim(merge('upper_divertor_inner', &
                'lower_divertor_outer', isupper)) // '_target) - ' // &
                'face region 4 will be empty'
        end if

        simgrid%face%reg = 0
        nfreg = 0
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) then
                tfc = GetFaceCell(simgrid%face, i)
                c = tfc(1)
                if (cvkind(c) == KCORE) then
                    ! Core boundary
                    simgrid%face%reg(i) = 8
                else
                    ! Nearest wall structure (trusted only when the
                    ! face actually sits on the wall)
                    xm = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%x(simgrid%face%vert(i, 2)))
                    ym = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%y(simgrid%face%vert(i, 2)))
                    k = 0
                    dmin = posinfval_R8()
                    do j = 1, int(vessel%nstructures, kind=I8)
                        d = StructureDistance(vessel, xm, ym, j)
                        if (d < dmin) then
                            dmin = d
                            k = j
                        end if
                    end do
                    if (dmin > dtolwall) k = 0
                    if (k > 0) then
                        if (istarget1(k) .and. istarget4(k)) then
                            ! Both targets share this wall structure
                            ! (e.g. both legs strike one flat wall
                            ! run): split it at the divertor (PFR)
                            ! cut, which the volume regions already
                            ! encode - West divertor cells (3) face
                            ! the West target, East cells (4) the
                            ! East one. Outside the divertor volumes
                            ! fall back to the side of the X-point ->
                            ! core axis. The PFR wall regions between
                            ! the targets then do not exist.
                            if (simgrid%cell%reg(c) == 3) then
                                simgrid%face%reg(i) = 1
                            elseif (simgrid%cell%reg(c) == 4) then
                                simgrid%face%reg(i) = 4
                            elseif ((ux*(ym - yxp) - &
                                uy*(xm - xxp)) > 0.0_R8) then
                                simgrid%face%reg(i) = 1
                            else
                                simgrid%face%reg(i) = 4
                            end if
                        elseif (istarget1(k)) then
                            simgrid%face%reg(i) = 1
                        elseif (istarget4(k)) then
                            simgrid%face%reg(i) = 4
                        end if
                    end if
                    if (simgrid%face%reg(i) == 0) then
                        if (cvkind(c) == KPFR) then
                            ! PFR wall, West/East half
                            if (simgrid%cell%reg(c) == 3) then
                                simgrid%face%reg(i) = 7
                            else
                                simgrid%face%reg(i) = 9
                            end if
                        else
                            ! Main-chamber wall: divertor SOL parts
                            ! keep their own regions
                            select case (simgrid%cell%reg(c))
                            case (3)
                                simgrid%face%reg(i) = 11
                            case (4)
                                simgrid%face%reg(i) = 13
                            case default
                                simgrid%face%reg(i) = 12
                            end select
                        end if
                    end if
                end if
            else
                ! Internal faces: separatrix, then the four cuts
                lbl = simgrid%face%TMfacelabel(i)
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) then
                    if (issepTM(lbl)) simgrid%face%reg(i) = 10
                end if
                if ((simgrid%face%reg(i) == 0) .and. &
                    (chainreg(i) /= 0)) then
                    simgrid%face%reg(i) = chainreg(i)
                end if
            end if
            nfreg(simgrid%face%reg(i)) = &
                nfreg(simgrid%face%reg(i)) + 1
        end do
        print *, 'SetSOLPSRegionsSN: fcReg counts (regions ' // &
            '1..13): ', nfreg(1:13)

        ! Flux tube regions
        !==================
        ! 1 = core, 2 = SOL (main + divertor), 3 = PFR
        nmix = 0
        nftreg = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            k = fd%fluxtubecells(fd%fluxtubecellsP(i, 1))
            fd%fluxtuberegID(i) = cvkind(k)
            nftreg(cvkind(k)) = nftreg(cvkind(k)) + 1
            do lbl = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (cvkind(fd%fluxtubecells(lbl)) /= cvkind(k)) &
                    nmix = nmix + 1
            end do
        end do
        end associate
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsSN: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL/PFR kinds ' // &
                '(tube crossing the separatrix?)'
        end if
        print *, 'SetSOLPSRegionsSN: ftReg: ', nftreg(1), &
            ' core tubes, ', nftreg(2), ' SOL tubes, ', nftreg(3), &
            ' PFR tubes'

    contains

        ! Flood fill of the cells of kind kindval, without crossing
        ! the chain faces tagged breg1/breg2, appending the new
        ! components to comp/ncomp (host-associated)
        subroutine FloodFill(kindval, breg1, breg2)
            integer(I8), intent(in)     :: kindval, breg1, breg2
            integer(I8)                 :: i0, jj, kk, cq, cn
            integer(I8), allocatable    :: cfl(:), tfcl(:)
            do i0 = 1, simgrid%cell%ntot
                if (cvkind(i0) /= kindval) cycle
                if (comp(i0) /= 0) cycle
                ncomp = ncomp + 1
                comp(i0) = ncomp
                nqh = 1
                nqt = 1
                queue(1) = i0
                do while (nqh <= nqt)
                    cq = queue(nqh)
                    nqh = nqh + 1
                    cfl = GetCellFace(simgrid%cell, cq)
                    do jj = 1, size(cfl)
                        if (simgrid%face%BF(cfl(jj))) cycle
                        if ((chainreg(cfl(jj)) == breg1) .or. &
                            (chainreg(cfl(jj)) == breg2)) cycle
                        tfcl = GetFaceCell(simgrid%face, cfl(jj))
                        do kk = 1, size(tfcl)
                            cn = tfcl(kk)
                            if (cn == cq) cycle
                            if (cvkind(cn) /= kindval) cycle
                            if (comp(cn) /= 0) cycle
                            comp(cn) = ncomp
                            nqt = nqt + 1
                            queue(nqt) = cn
                        end do
                    end do
                end do
            end do
        end subroutine

    end subroutine


    ! SOLPS region conventions for the DN (double null) catalog
    ! topologies: CDN, DDN bottom and DDN top
    subroutine SetSOLPSRegionsDN(simgrid, topomesh, vessel, isddn)

        ! Description
        !============
        ! Assign the volume and face regions that SOLPS expects for a
        ! double-null case (NREG = 8; FREG = 26 for CDN, 27 for DDN),
        ! 
        ! on the quadrilateral-equivalent picture of the unstructured
        ! Goat grid. The divertor numbering goes clockwise starting
        ! from the lower inner target. Only the innermost X-point per
        ! lower/upper side takes part in the construction (consistent
        ! with ClassifyBasicSOLPSCatalogTopology).
        !
        ! Each of the two X-points contributes 4 separatrix branches
        ! and 4 radial chains of non-aligned grid faces walked
        ! outwards from the X-point vertex: into the core-directed
        ! sector (core cuts 9/11; for a DDN the secondary X-point's
        ! chain first crosses the SOL band between the two
        ! separatrices - those faces are region 13, the connection
        ! between right and left SOL - and continues through the
        ! core), into the PFR (PFR cuts 12/10), and into the two SOL
        ! sectors (divertor entrances 2/7 below, 3/6 above).
        !
        ! Volume regions (cvReg):
        !   1 = left (inboard) core        5 = right (outboard) core
        !   2 = left SOL                   6 = right SOL
        !   3 = lower inner divertor       4 = upper inner divertor
        !   7 = upper outer divertor       8 = lower outer divertor
        ! (divertor = divertor SOL + PFR half). The cores are the
        ! flood-fill components of the core cells with the two core
        ! cuts as barriers; everything else comes from a single
        ! flood fill of the non-core cells with all chains as
        ! barriers: the components adjacent to the core across the
        ! separatrix are the left/right SOL, the remaining ones the
        ! four divertors (sorted by side and radius).
        !
        ! Face regions (fcReg), FREG = 26 (CDN) / 27 (DDN); the
        ! DDN-only region 13 shifts every region >= 13 up by one:
        !   1/4/5/8   = lower inner / upper inner / upper outer /
        !               lower outer target (wall faces on the
        !               goat.vessel.solps_*_divertor_*_target
        !               structures). The inner and outer list of a
        !               side may share a structure (both legs
        !               striking one flat wall run): its faces are
        !               then split between the two target regions at
        !               that side's PFR cut (via the volume regions),
        !               and the PFR wall regions in between do not
        !               exist.
        !   2/3/6/7   = entrances to the lower inner / upper inner /
        !               upper outer / lower outer divertor
        !   9/11      = lower/upper core cut
        !   10/12     = upper/lower PFR cut
        !   13        = (DDN only) connection between right and left
        !               SOL: the between-separatrix segment of the
        !               secondary X-point's core-directed chain
        !   13/15/20/22 (+1 for DDN) = lower inner / upper inner /
        !               upper outer / lower outer PFR wall boundary
        !   14/21 (+1) = left/right core boundary
        !   16/23 (+1) = left/right separatrix (the core-enclosing
        !               separatrix only, split by the X-points/cuts)
        !   17/19/24/26 (+1) = lower inner / upper inner / upper
        !               outer / lower outer divertor main wall
        !   18/25 (+1) = left/right SOL main wall boundary
        !
        ! Flux tube regions (ftReg): 1 = core, 2 = SOL, 3 = PFR
        ! (PFR = divertor cells on the core-side flux of their own
        ! X-point, decided per zone).
        !
        ! Empty target inputs are tolerated (warning): the target
        ! face region stays empty and its wall faces fall to the
        ! PFR/divertor/SOL wall regions. Faces are attributed to a
        ! wall structure only within dtolwall.

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel
        logical, intent(in)             :: isddn

        ! Parameters: cell kinds and wall attribution tolerance
        integer(I8), parameter          :: KCORE = 1, KSOL = 2, &
            KPFR = 3
        real(R8), parameter             :: dtolwall = 5.0e-3_R8

        ! Auxiliary
        integer(I8), allocatable        :: xps(:), op(:), cc(:), &
            ccf(:), tmpi(:), tcv(:), tfc(:), vf(:), cf(:), zn(:), &
            zoneid(:), cvkind(:), comp(:), queue(:), chainreg(:), &
            compreg(:), compn(:), compentrance(:)
        logical, allocatable            :: iscorezone(:), &
            issepTM(:), iscoreadj(:), istLI(:), istLO(:), &
            istUI(:), istUO(:)
        real(R8), allocatable           :: zpsi(:), compx(:), &
            compy(:)
        integer(I8)                     :: i, j, k, ich, iv, w, f, &
            fprev, lbl, c, c2, nxp, ncomp, ncorecomp, nqh, nqt, &
            ilow, iup, off, nmix, n13c, kside, ivX(2), fstart(2, 4), &
            chlen(2, 4), termvx(2, 4), ncand(2, 4), nfreg(0:27), &
            ncvreg(8), nftreg(3), nprob1
        real(R8)                        :: xc, yc, fcore, psicore, &
            s, d, dmin, score, dirx, diry, dp, xm, ym, vx, vy, zm, &
            xxp(2), yxp(2), psiX(2)
        character(14)                   :: chname(4)
        character(5)                    :: sidename(2)

        chname = ['core-directed ', 'PFR cut       ', &
            'inner entrance', 'outer entrance']
        sidename = ['lower', 'upper']
        off = merge(1_I8, 0_I8, isddn)
        print *, 'SetSOLPSRegionsDN: applying SOLPS double-null ' // &
            'regions (' // trim(merge('DDN', 'CDN', isddn)) // ')'

        ! Reference center and core-side flux value
        !==========================================
        ! From the O-point vertex when the axis is meshed, otherwise
        ! from the core-cell centroid and the core boundary faces
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            print *, 'SetSOLPSRegionsDN: no core topomesh cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        op = topomesh%GetOPointIDs()
        if (size(op) == 1) then
            xc = topomesh%vert%x(op(1))
            yc = topomesh%vert%y(op(1))
            fcore = topomesh%vert%fval(op(1))
        else
            allocate(tmpi(0))
            do i = 1, size(cc)
                tmpi = [tmpi, topomesh%cell%GetVert(cc(i))]
            end do
            xc = sum(topomesh%vert%x(tmpi))/real(size(tmpi), kind=R8)
            yc = sum(topomesh%vert%y(tmpi))/real(size(tmpi), kind=R8)
            deallocate(tmpi)
            ccf = topomesh%GetCoreFaceIDs()
            if (size(ccf) == 0) then
                print *, 'SetSOLPSRegionsDN: no core-side flux ' // &
                    'reference - leaving the generic regions in place'
                return
            end if
            fcore = 0.0_R8
            do i = 1, size(ccf)
                fcore = fcore + &
                    0.5_R8*(topomesh%vert%fval( &
                            topomesh%face%vert(ccf(i), 1)) + &
                            topomesh%vert%fval( &
                            topomesh%face%vert(ccf(i), 2)))
            end do
            fcore = fcore/real(size(ccf), kind=R8)
        end if

        ! The two X-points: innermost per side
        !=====================================
        xps = topomesh%GetXPointIDs()
        ! Narrow-grid rule: an X-point whose flux surface was removed
        ! from the delivered grid by the narrow-grid target filter
        ! (orphaned grid vertex, fieldlineID 0) does not exist for the
        ! SOLPS region conventions - drop it from the working list. On
        ! a complete (wide) grid every X-point stays.
        if (size(xps) > 0) then
            if (any(simgrid%vert%fieldlineID(xps) <= 0)) then
                print *, 'SetSOLPSRegionsDN: ignoring ', &
                    count(simgrid%vert%fieldlineID(xps) <= 0), &
                    ' X-point(s) without a flux surface in the ' // &
                    'filtered grid'
                xps = pack(xps, simgrid%vert%fieldlineID(xps) > 0)
            end if
        end if
        nxp = size(xps)
        ilow = 0
        iup = 0
        vx = 0.0_R8
        vy = 0.0_R8
        do i = 1, nxp
            d = abs(topomesh%vert%fval(xps(i)) - fcore)
            if (topomesh%vert%y(xps(i)) >= yc) then
                if ((iup == 0) .or. (d < vy)) then
                    iup = i
                    vy = d
                end if
            else
                if ((ilow == 0) .or. (d < vx)) then
                    ilow = i
                    vx = d
                end if
            end if
        end do
        if ((ilow == 0) .or. (iup == 0)) then
            print *, 'SetSOLPSRegionsDN: could not find one ' // &
                'X-point per side - leaving the generic regions ' // &
                'in place'
            return
        end if
        xxp(1) = topomesh%vert%x(xps(ilow))
        yxp(1) = topomesh%vert%y(xps(ilow))
        xxp(2) = topomesh%vert%x(xps(iup))
        yxp(2) = topomesh%vert%y(xps(iup))

        ! Fail-fast: each primary X-point's strike legs must land on its
        ! side's declared divertor targets (all four target lists must
        ! be given). The lower primary feeds the lower inner/outer, the
        ! upper primary the upper inner/outer.
        ! Check every X-point (primaries + any snowflake companions);
        ! the guard keeps only the real divertor legs (two per X-point).
        call GuardPrimaryStrikes(topomesh, vessel, dtolwall, &
            'SetSOLPSRegionsDN', xps, &
            [1_I8, 2_I8, 3_I8, 4_I8], nprob1)
        if (nprob1 > 0) then
            call gdErrorHandler('SetSOLPSRegionsDN: primary strike ' // &
                'points are not all on declared divertor targets ' // &
                '(see the messages above) - set/split the ' // &
                'goat.vessel.solps_{lower,upper}_divertor_' // &
                '{inner,outer}_target lists so every primary strike ' // &
                'lands on target wall')
        end if

        ! Grid vertices sitting on the X-points
        do k = 1, 2
            ivX(k) = 1
            dmin = posinfval_R8()
            do i = 1, simgrid%vert%ntot
                d = (simgrid%vert%x(i) - xxp(k))**2 + &
                    (simgrid%vert%y(i) - yxp(k))**2
                if (d < dmin) then
                    dmin = d
                    ivX(k) = i
                end if
            end do
            if (sqrt(dmin) > 1.0e-6_R8) then
                print *, 'SetSOLPSRegionsDN: WARNING: nearest ' // &
                    'grid vertex to the ' // sidename(k) // &
                    ' X-point is ', sqrt(dmin), ' m away'
            end if
            psiX(k) = simgrid%vert%psi(ivX(k))
        end do

        ! Zones and core cells
        !=====================
        ! Snapshot of the pre-overwrite regions (the topomesh cells,
        ! which never straddle a separatrix): zone mean flux for the
        ! PFR tests, core zones from the region marks
        allocate(zoneid(simgrid%cell%ntot))
        zoneid = simgrid%cell%reg
        allocate(zpsi(maxval(zoneid)), zn(maxval(zoneid)), &
            iscorezone(maxval(zoneid)))
        zpsi = 0.0_R8
        zn = 0
        do i = 1, simgrid%cell%ntot
            zpsi(zoneid(i)) = zpsi(zoneid(i)) + simgrid%cell%psi(i)
            zn(zoneid(i)) = zn(zoneid(i)) + 1
        end do
        iscorezone = .false.
        psicore = 0.0_R8
        k = 0
        do j = 1, size(zn)
            if (zn(j) == 0) cycle
            if (mod(j - SOLPScoreregID, SOLPScoreregIDincr) == 0) &
                then
                iscorezone(j) = .true.
                psicore = psicore + zpsi(j)
                k = k + zn(j)
            end if
        end do
        if (k == 0) then
            print *, 'SetSOLPSRegionsDN: no core grid cells - ' // &
                'leaving the generic regions in place'
            return
        end if
        s = sign(1.0_R8, psiX(1) - psicore/real(k, kind=R8))
        allocate(cvkind(simgrid%cell%ntot))
        do i = 1, simgrid%cell%ntot
            cvkind(i) = merge(KCORE, KSOL, iscorezone(zoneid(i)))
        end do

        ! The four radial rays leaving each X-point vertex
        !=================================================
        ! Classified by the sector they walk into: core-side flux
        ! rays split into core-directed (towards the core centroid)
        ! and PFR (away from it); SOL-side rays split into inner and
        ! outer by the side of the X-point -> core direction
        fstart = 0
        ncand = 0
        do k = 1, 2
            vx = xc - xxp(k)
            vy = yc - yxp(k)
            vf = GetVertFace(simgrid%vert, ivX(k))
            do j = 1, size(vf)
                f = vf(j)
                if (simgrid%face%BF(f)) cycle
                if (simgrid%face%aligned(f) == 1) cycle
                if (simgrid%face%vert(f, 1) == ivX(k)) then
                    w = simgrid%face%vert(f, 2)
                else
                    w = simgrid%face%vert(f, 1)
                end if
                dirx = simgrid%vert%x(w) - xxp(k)
                diry = simgrid%vert%y(w) - yxp(k)
                dp = s*(simgrid%vert%psi(w) - psiX(k))
                if (dp < 0.0_R8) then
                    if ((dirx*vx + diry*vy) > 0.0_R8) then
                        ich = 1
                    else
                        ich = 2
                    end if
                else
                    ! Inner = left of u for the lower X-point,
                    ! right of u for the upper one
                    if (((vx*diry - vy*dirx) > 0.0_R8) .eqv. &
                        (k == 1)) then
                        ich = 3
                    else
                        ich = 4
                    end if
                end if
                ncand(k, ich) = ncand(k, ich) + 1
                if (fstart(k, ich) == 0) fstart(k, ich) = f
            end do
            print *, 'SetSOLPSRegionsDN: ' // sidename(k) // &
                ' X-point grid vertex ', ivX(k), ' at ', xxp(k), &
                yxp(k), ': radial ray candidates (core-directed/' // &
                'PFR/inner/outer): ', ncand(k, :)
            do ich = 1, 4
                if (fstart(k, ich) == 0) then
                    print *, 'SetSOLPSRegionsDN: WARNING: no ' // &
                        trim(chname(ich)) // ' ray at the ' // &
                        sidename(k) // ' X-point'
                end if
            end do
        end do

        ! Walk the eight radial chains
        !=============================
        ! Face tags: entrances 2/7 (lower inner/outer), 3/6 (upper
        ! inner/outer); PFR cuts 12 (lower), 10 (upper); the
        ! core-directed chains are tagged per face: 9/11 (lower/
        ! upper core cut) next to core cells, 13 (between the two
        ! separatrices, DDN only) next to SOL cells
        allocate(chainreg(simgrid%face%ntot))
        chainreg = 0
        chlen = 0
        termvx = 0
        do k = 1, 2
            do ich = 1, 4
                if (fstart(k, ich) == 0) cycle
                iv = ivX(k)
                f = fstart(k, ich)
                do
                    select case (ich)
                    case (1)
                        tfc = GetFaceCell(simgrid%face, f)
                        if (any(iscorezone(zoneid(tfc)))) then
                            chainreg(f) = merge(9_I8, 11_I8, k == 1)
                        else
                            chainreg(f) = 13
                        end if
                    case (2)
                        chainreg(f) = merge(12_I8, 10_I8, k == 1)
                    case (3)
                        chainreg(f) = merge(2_I8, 3_I8, k == 1)
                    case (4)
                        chainreg(f) = merge(7_I8, 6_I8, k == 1)
                    end select
                    chlen(k, ich) = chlen(k, ich) + 1
                    if (chlen(k, ich) > simgrid%face%ntot) then
                        print *, 'SetSOLPSRegionsDN: WARNING: ' // &
                            'runaway chain walk - aborted'
                        exit
                    end if
                    if (simgrid%face%vert(f, 1) == iv) then
                        w = simgrid%face%vert(f, 2)
                    else
                        w = simgrid%face%vert(f, 1)
                    end if
                    if (simgrid%vert%BV(w)) then
                        termvx(k, ich) = w
                        exit
                    end if
                    vx = simgrid%vert%x(w) - simgrid%vert%x(iv)
                    vy = simgrid%vert%y(w) - simgrid%vert%y(iv)
                    fprev = f
                    f = 0
                    score = -posinfval_R8()
                    cf = GetVertFace(simgrid%vert, w)
                    do j = 1, size(cf)
                        if (cf(j) == fprev) cycle
                        if (simgrid%face%BF(cf(j))) cycle
                        if (simgrid%face%aligned(cf(j)) == 1) cycle
                        if (chainreg(cf(j)) /= 0) cycle
                        if (simgrid%face%vert(cf(j), 1) == w) then
                            c2 = simgrid%face%vert(cf(j), 2)
                        else
                            c2 = simgrid%face%vert(cf(j), 1)
                        end if
                        dirx = simgrid%vert%x(c2) - simgrid%vert%x(w)
                        diry = simgrid%vert%y(c2) - simgrid%vert%y(w)
                        d = (dirx*vx + diry*vy)/ &
                            max(sqrt(dirx**2 + diry**2), &
                            tiny(1.0_R8))
                        if (d > score) then
                            score = d
                            f = cf(j)
                        end if
                    end do
                    if (f == 0) then
                        print *, 'SetSOLPSRegionsDN: WARNING: ' // &
                            trim(chname(ich)) // ' chain of the ' // &
                            sidename(k) // ' X-point stopped at ' // &
                            'internal vertex ', w, ' after ', &
                            chlen(k, ich), ' faces'
                        termvx(k, ich) = w
                        exit
                    end if
                    iv = w
                end do
                if (termvx(k, ich) /= 0) then
                    print *, 'SetSOLPSRegionsDN: ' // sidename(k) // &
                        ' ' // trim(chname(ich)) // ' chain: ', &
                        chlen(k, ich), ' faces, ends at ', &
                        simgrid%vert%x(termvx(k, ich)), &
                        simgrid%vert%y(termvx(k, ich))
                end if
            end do
        end do

        ! Volume regions
        !===============
        ! Flood fill of the core cells with the core cuts as
        ! barriers (left/right core), then of all the non-core
        ! cells with all the chains as barriers (left/right SOL +
        ! the four divertors)
        allocate(comp(simgrid%cell%ntot), queue(simgrid%cell%ntot))
        comp = 0
        ncomp = 0
        call FloodFill(.true.)
        ncorecomp = ncomp
        call FloodFill(.false.)
        if (ncorecomp /= 2) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', ncorecomp, &
                ' core components instead of 2 (incomplete core ' // &
                'cuts?)'
        end if
        if ((ncomp - ncorecomp) /= 6) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', &
                ncomp - ncorecomp, ' non-core components instead ' // &
                'of 6 (incomplete chains?)'
        end if

        ! Component classification: centroids, core adjacency
        allocate(compx(ncomp), compy(ncomp), compn(ncomp), &
            compreg(ncomp), iscoreadj(ncomp), compentrance(ncomp))
        compx = 0.0_R8
        compy = 0.0_R8
        compn = 0
        do i = 1, simgrid%cell%ntot
            compx(comp(i)) = compx(comp(i)) + simgrid%cell%x(i)
            compy(comp(i)) = compy(comp(i)) + simgrid%cell%y(i)
            compn(comp(i)) = compn(comp(i)) + 1
        end do
        iscoreadj = .false.
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) cycle
            if (chainreg(i) /= 0) cycle
            tfc = GetFaceCell(simgrid%face, i)
            if (size(tfc) /= 2) cycle
            if ((cvkind(tfc(1)) == KCORE) .neqv. &
                (cvkind(tfc(2)) == KCORE)) then
                if (cvkind(tfc(1)) == KCORE) then
                    iscoreadj(comp(tfc(2))) = .true.
                else
                    iscoreadj(comp(tfc(1))) = .true.
                end if
            end if
        end do

        ! Which entrance chain each component touches (on its divertor,
        ! non-core-adjacent side) -> the divertor it belongs to. This is what
        ! actually splits inner from outer (the SAME radial cut that makes the
        ! entrance/PFR face regions), robust to the two halves lying on the same
        ! radial side of the X-point (a centroid-radius test then misclassifies
        ! both as inner). Entrance chainreg: 2=lower inner, 3=upper inner,
        ! 6=upper outer, 7=lower outer -> volume regions 3/4/7/8 respectively.
        compentrance = 0
        do i = 1, simgrid%face%ntot
            k = chainreg(i)
            if ((k /= 2) .and. (k /= 3) .and. (k /= 6) .and. (k /= 7)) cycle
            tfc = GetFaceCell(simgrid%face, i)
            do j = 1, size(tfc)
                if (cvkind(tfc(j)) == KCORE) cycle
                if (iscoreadj(comp(tfc(j)))) cycle   ! the main-SOL side
                compentrance(comp(tfc(j))) = k       ! the divertor side
            end do
        end do

        ! Region of each component: cores 1/5 and SOLs 2/6 by the
        ! side of the reference center; divertors 3/4/7/8 by which entrance
        ! chain the component touches (the radial cut that also draws the
        ! entrance/PFR face regions -- this is what splits inner from outer).
        ! Only when a divertor component touches no entrance chain (should not
        ! happen for a well-formed grid) fall back to the side (below/above the
        ! center) and centroid radius vs the X-point.
        do c = 1, ncomp
            if (compn(c) == 0) then
                compreg(c) = 0
                cycle
            end if
            xm = compx(c)/real(compn(c), kind=R8)
            ym = compy(c)/real(compn(c), kind=R8)
            if (c <= ncorecomp) then
                compreg(c) = merge(1_I8, 5_I8, xm < xc)
            elseif (iscoreadj(c)) then
                compreg(c) = merge(2_I8, 6_I8, xm < xc)
            elseif (compentrance(c) == 2) then
                compreg(c) = 3_I8          ! lower inner divertor
            elseif (compentrance(c) == 3) then
                compreg(c) = 4_I8          ! upper inner divertor
            elseif (compentrance(c) == 6) then
                compreg(c) = 7_I8          ! upper outer divertor
            elseif (compentrance(c) == 7) then
                compreg(c) = 8_I8          ! lower outer divertor
            else
                kside = merge(1_I8, 2_I8, ym < yc)
                if (kside == 1) then
                    compreg(c) = merge(3_I8, 8_I8, xm < xxp(1))
                else
                    compreg(c) = merge(4_I8, 7_I8, xm < xxp(2))
                end if
            end if
        end do
        do i = 1, simgrid%cell%ntot
            simgrid%cell%reg(i) = compreg(comp(i))
        end do
        ncvreg = 0
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k >= 1) .and. (k <= 8)) ncvreg(k) = ncvreg(k) + 1
        end do
        print *, 'SetSOLPSRegionsDN: cvReg counts: left core ', &
            ncvreg(1), ', left SOL ', ncvreg(2), ', lower inner ' // &
            'div ', ncvreg(3), ', upper inner div ', ncvreg(4), &
            ', right core ', ncvreg(5), ', right SOL ', ncvreg(6), &
            ', upper outer div ', ncvreg(7), ', lower outer div ', &
            ncvreg(8)

        ! Cell kinds for the flux tube regions and the PFR walls:
        ! divertor cells on the core-side flux of their own X-point
        ! (per zone) are PFR
        do i = 1, simgrid%cell%ntot
            k = simgrid%cell%reg(i)
            if ((k == 3) .or. (k == 8)) then
                kside = 1
            elseif ((k == 4) .or. (k == 7)) then
                kside = 2
            else
                cycle
            end if
            zm = zpsi(zoneid(i))/real(zn(zoneid(i)), kind=R8)
            if (s*(zm - psiX(kside)) < 0.0_R8) cvkind(i) = KPFR
        end do

        ! Face regions
        !=============
        ! Separatrix topomesh faces: only the core-enclosing part
        ! (interfaces between a core and a non-core topomesh cell)
        allocate(issepTM(topomesh%face%ntot))
        issepTM = .false.
        do i = 1, topomesh%face%ntot
            if (topomesh%face%type(i) /= TMfacesepID) cycle
            tfc = topomesh%face%GetCell(i)
            if (size(tfc) == 2) then
                if (any(tfc(1) == cc) .neqv. any(tfc(2) == cc)) then
                    issepTM(i) = .true.
                end if
            end if
        end do

        ! Target structures from the Goat inputs
        allocate(istLI(vessel%nstructures), &
            istLO(vessel%nstructures), istUI(vessel%nstructures), &
            istUO(vessel%nstructures))
        istLI = .false.
        istLO = .false.
        istUI = .false.
        istUO = .false.
        do i = 1, int(vessel%nsolpsldi, kind=I8)
            istLI(vessel%solpsldiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsldo, kind=I8)
            istLO(vessel%solpsldoind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudi, kind=I8)
            istUI(vessel%solpsudiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudo, kind=I8)
            istUO(vessel%solpsudoind(i)) = .true.
        end do
        if (count(istLI) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_lower_divertor_inner_' // &
            'target empty - face region 1 will be empty'
        if (count(istUI) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_upper_divertor_inner_' // &
            'target empty - face region 4 will be empty'
        if (count(istUO) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_upper_divertor_outer_' // &
            'target empty - face region 5 will be empty'
        if (count(istLO) == 0) print *, 'SetSOLPSRegionsDN: ' // &
            'WARNING: goat.vessel.solps_lower_divertor_outer_' // &
            'target empty - face region 8 will be empty'

        simgrid%face%reg = 0
        nfreg = 0
        n13c = 0
        do i = 1, simgrid%face%ntot
            if (simgrid%face%BF(i)) then
                tfc = GetFaceCell(simgrid%face, i)
                c = tfc(1)
                if (cvkind(c) == KCORE) then
                    ! Left/right core boundary
                    simgrid%face%reg(i) = &
                        merge(14_I8, 21_I8, &
                        simgrid%cell%reg(c) == 1) + off
                else
                    ! Nearest wall structure (trusted only when the
                    ! face actually sits on the wall)
                    xm = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%x(simgrid%face%vert(i, 2)))
                    ym = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) &
                        + simgrid%vert%y(simgrid%face%vert(i, 2)))
                    k = 0
                    dmin = posinfval_R8()
                    do j = 1, int(vessel%nstructures, kind=I8)
                        d = StructureDistance(vessel, xm, ym, j)
                        if (d < dmin) then
                            dmin = d
                            k = j
                        end if
                    end do
                    if (dmin > dtolwall) k = 0
                    if (k > 0) then
                        if (istLI(k) .and. istLO(k)) then
                            ! Both lower targets share this wall
                            ! structure: split it at the lower PFR
                            ! cut, which the volume regions already
                            ! encode (3 = lower inner divertor,
                            ! 8 = lower outer; left/right SOL fall
                            ! to the matching side). The PFR wall
                            ! regions between the targets then do
                            ! not exist.
                            select case (simgrid%cell%reg(c))
                            case (3)
                                simgrid%face%reg(i) = 1
                            case (8)
                                simgrid%face%reg(i) = 8
                            case (6)
                                simgrid%face%reg(i) = 8
                            case default
                                simgrid%face%reg(i) = 1
                            end select
                        elseif (istUI(k) .and. istUO(k)) then
                            ! Both upper targets share this wall
                            ! structure: split it at the upper PFR
                            ! cut (4 = upper inner divertor,
                            ! 7 = upper outer)
                            select case (simgrid%cell%reg(c))
                            case (4)
                                simgrid%face%reg(i) = 4
                            case (7)
                                simgrid%face%reg(i) = 5
                            case (6)
                                simgrid%face%reg(i) = 5
                            case default
                                simgrid%face%reg(i) = 4
                            end select
                        else
                            if (istLI(k)) simgrid%face%reg(i) = 1
                            if (istUI(k)) simgrid%face%reg(i) = 4
                            if (istUO(k)) simgrid%face%reg(i) = 5
                            if (istLO(k)) simgrid%face%reg(i) = 8
                        end if
                    end if
                    if (simgrid%face%reg(i) == 0) then
                        select case (simgrid%cell%reg(c))
                        case (2)
                            simgrid%face%reg(i) = 18 + off
                        case (6)
                            simgrid%face%reg(i) = 25 + off
                        case (3)
                            simgrid%face%reg(i) = merge(13_I8, &
                                17_I8, cvkind(c) == KPFR) + off
                        case (4)
                            simgrid%face%reg(i) = merge(15_I8, &
                                19_I8, cvkind(c) == KPFR) + off
                        case (7)
                            simgrid%face%reg(i) = merge(20_I8, &
                                24_I8, cvkind(c) == KPFR) + off
                        case (8)
                            simgrid%face%reg(i) = merge(22_I8, &
                                26_I8, cvkind(c) == KPFR) + off
                        case default
                            simgrid%face%reg(i) = 18 + off
                        end select
                    end if
                end if
            else
                ! Internal faces: left/right separatrix by the core
                ! component behind them, then the chains
                lbl = simgrid%face%TMfacelabel(i)
                if ((lbl >= 1) .and. (lbl <= topomesh%face%ntot)) &
                    then
                    if (issepTM(lbl)) then
                        tfc = GetFaceCell(simgrid%face, i)
                        c = tfc(1)
                        if (cvkind(c) /= KCORE) then
                            if (size(tfc) == 2) c = tfc(2)
                        end if
                        simgrid%face%reg(i) = merge(16_I8, 23_I8, &
                            simgrid%cell%reg(c) == 1) + off
                    end if
                end if
                if ((simgrid%face%reg(i) == 0) .and. &
                    (chainreg(i) /= 0)) then
                    if ((chainreg(i) == 13) .and. .not. isddn) then
                        ! No between-separatrix connection exists in
                        ! a CDN: leave unassigned and warn
                        n13c = n13c + 1
                    else
                        simgrid%face%reg(i) = chainreg(i)
                    end if
                end if
            end if
            nfreg(simgrid%face%reg(i)) = &
                nfreg(simgrid%face%reg(i)) + 1
        end do
        if (n13c > 0) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', n13c, &
                ' between-separatrix chain faces in a CDN ' // &
                '(imperfectly connected null?) - left unassigned'
        end if
        print *, 'SetSOLPSRegionsDN: fcReg counts (regions 1..', &
            26 + off, '): ', nfreg(1:26 + off)

        ! Flux tube regions
        !==================
        ! 1 = core, 2 = SOL (main + divertor), 3 = PFR
        nmix = 0
        nftreg = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            k = fd%fluxtubecells(fd%fluxtubecellsP(i, 1))
            fd%fluxtuberegID(i) = cvkind(k)
            nftreg(cvkind(k)) = nftreg(cvkind(k)) + 1
            do lbl = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (cvkind(fd%fluxtubecells(lbl)) /= cvkind(k)) &
                    nmix = nmix + 1
            end do
        end do
        end associate
        if (nmix > 0) then
            print *, 'SetSOLPSRegionsDN: WARNING: ', nmix, &
                ' flux tube cells with mixed core/SOL/PFR kinds ' // &
                '(tube crossing a separatrix?)'
        end if
        print *, 'SetSOLPSRegionsDN: ftReg: ', nftreg(1), &
            ' core tubes, ', nftreg(2), ' SOL tubes, ', nftreg(3), &
            ' PFR tubes'

    contains

        ! Flood fill of the core (docore) or non-core (.not.docore)
        ! cells, without crossing any chain face, appending the new
        ! components to comp/ncomp (host-associated)
        subroutine FloodFill(docore)
            logical, intent(in)         :: docore
            integer(I8)                 :: i0, jj, kk, cq, cn
            integer(I8), allocatable    :: cfl(:), tfcl(:)
            do i0 = 1, simgrid%cell%ntot
                if ((cvkind(i0) == KCORE) .neqv. docore) cycle
                if (comp(i0) /= 0) cycle
                ncomp = ncomp + 1
                comp(i0) = ncomp
                nqh = 1
                nqt = 1
                queue(1) = i0
                do while (nqh <= nqt)
                    cq = queue(nqh)
                    nqh = nqh + 1
                    cfl = GetCellFace(simgrid%cell, cq)
                    do jj = 1, size(cfl)
                        if (simgrid%face%BF(cfl(jj))) cycle
                        if (chainreg(cfl(jj)) /= 0) cycle
                        tfcl = GetFaceCell(simgrid%face, cfl(jj))
                        do kk = 1, size(tfcl)
                            cn = tfcl(kk)
                            if (cn == cq) cycle
                            if ((cvkind(cn) == KCORE) .neqv. &
                                docore) cycle
                            if (comp(cn) /= 0) cycle
                            comp(cn) = ncomp
                            nqt = nqt + 1
                            queue(nqt) = cn
                        end do
                    end do
                end do
            end do
        end subroutine

    end subroutine

    ! Fail-fast guard: every strike leg of the given primary X-point(s)
    ! must land (within tol of the wall) on a structure declared as a
    ! divertor target in one of the required target lists (integer codes
    ! 1 = lower inner, 2 = lower outer, 3 = upper inner, 4 = upper
    ! outer). ALL the required lists must be non-empty. The allowed set
    ! is the UNION of the required lists, NOT per-side: a disconnected
    ! double null's secondary separatrix strikes BOTH divertors, so an
    ! upper X-point's legs can legitimately land on the lower targets
    ! (this matches the plan's "one of the four target lists" rule). Two
    ! lists may point to the same single wall structure (both legs on
    ! one flat run), but every required list must be given. This forces
    ! the user to declare the target segments in preprocessing so the
    ! target regions are attributed correctly - a deliberate exception
    ! to the region writers' otherwise no-fail philosophy (see the
    ! SECONDARY_XPOINT_REGIONS_PLAN, preliminary task). nproblems returns
    ! the number of problems (missing lists + off-target strikes); the
    ! caller raises a single error if > 0, so the user sees every
    ! offending strike (all divertors) in one run.
    subroutine GuardPrimaryStrikes(topomesh, vessel, tol, writer, &
        xpverts, codes, nproblems)

        ! Arguments
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel
        real(R8), intent(in)            :: tol
        character(*), intent(in)        :: writer
        integer(I8), intent(in)         :: xpverts(:), codes(:)
        integer(I8), intent(out)        :: nproblems

        ! Auxiliary
        integer(I8), allocatable        :: sp(:), spx(:), candsp(:), &
            usedsp(:)
        real(R8), allocatable           :: candang(:)
        logical, allocatable            :: allowed(:)
        integer(I8)                     :: i, k, ic, kbest, nstrike, &
            nviol, ncand, nused, imin, imax
        real(R8)                        :: xs, ys, d, dmin, xX, yX, &
            xcen, ycen, ux, uy
        logical                         :: lists_ok, ontarget, okc

        nproblems = 0

        ! Allowed target structures = union of the required lists (built
        ! with count-guarded loops, so an unallocated/empty list is safe).
        ! Every required list must be non-empty.
        allocate(allowed(vessel%nstructures))
        allowed = .false.
        lists_ok = .true.
        do ic = 1, size(codes)
            call MarkList(codes(ic))
            if (ListCount(codes(ic)) <= 0) lists_ok = .false.
        end do

        ! Strike legs of the given primary X-point(s). When anything is
        ! wrong, report each strike with its nearest structure so the
        ! user knows which structures to declare as targets.
        !
        ! Only strikes on divertor legs FORMED BY their own X-point are
        ! used for this target-attribution criterion. GOAT's strike
        ! algorithm (GetStrikePointIDs) is unchanged and still returns
        ! every wall-touching separatrix boundary point; but a secondary
        ! separatrix can graze the main wall on its core-side arms and
        ! re-enter, producing extra "strike" points that are NOT divertor
        ! legs of their X-point (e.g. an HFS-midplane touch). Such a leg
        ! is recognized geometrically: a genuine divertor leg leaves the
        ! X-point towards the divertor, i.e. the strike lies on the far
        ! side of the X-point from the core center; the grazing arms point
        ! back towards the core and are dropped. This leaves the two real
        ! divertor legs per X-point.
        sp = topomesh%GetStrikePointIDs()
        spx = topomesh%GetStrikePointXPointIDs()
        call SOLPSCoreCenter(topomesh, xcen, ycen, okc)
        allocate(candsp(size(sp)), candang(size(sp)), usedsp(size(sp)))
        nstrike = 0
        nviol = 0

        ! One X-point at a time: gather its divertor-side strikes and keep
        ! the two divertor legs. NOTE (edge cases still to handle): an SF+
        ! X-point sitting in a PFR contributes no real target strikes, and
        ! an SF0 X-point exactly on the original leg keeps its central
        ! strike as the original one - neither is treated specially yet.
        do ic = 1, size(xpverts)
            xX = topomesh%vert%x(xpverts(ic))
            yX = topomesh%vert%y(xpverts(ic))
            ux = xX - xcen          ! divertor direction (away from core)
            uy = yX - ycen
            ncand = 0
            do i = 1, size(sp)
                if (spx(i) /= xpverts(ic)) cycle
                xs = topomesh%vert%x(sp(i))
                ys = topomesh%vert%y(sp(i))
                ! drop core-side arms (grazing touches that point back
                ! towards the core, not into the divertor)
                if (okc) then
                    if (((xs - xX)*ux + (ys - yX)*uy) <= 0.0_R8) cycle
                end if
                ncand = ncand + 1
                candsp(ncand) = sp(i)
                if (okc) then
                    candang(ncand) = atan2( &
                        ux*(ys - yX) - uy*(xs - xX), &
                        ux*(xs - xX) + uy*(ys - yX))
                else
                    candang(ncand) = 0.0_R8
                end if
            end do

            ! Keep only the two divertor legs = the two angular extremes
            ! (the legs bounding the divertor/new PFR). Any central leg
            ! between them is the separatrix continuation inside the new
            ! PFR (a snowflake), not a target strike, and is dropped.
            if (ncand > 2) then
                imin = 1
                imax = 1
                do i = 2, ncand
                    if (candang(i) < candang(imin)) imin = i
                    if (candang(i) > candang(imax)) imax = i
                end do
                nused = 2
                usedsp(1) = candsp(imin)
                usedsp(2) = candsp(imax)
            else
                nused = ncand
                do i = 1, ncand
                    usedsp(i) = candsp(i)
                end do
            end if

            ! Validate the two divertor legs against the declared targets
            do i = 1, nused
                xs = topomesh%vert%x(usedsp(i))
                ys = topomesh%vert%y(usedsp(i))
                nstrike = nstrike + 1
                kbest = 0
                dmin = posinfval_R8()
                do k = 1, int(vessel%nstructures, kind=I8)
                    d = StructureDistance(vessel, xs, ys, k)
                    if (d < dmin) then
                        dmin = d
                        kbest = k
                    end if
                end do
                ontarget = (kbest > 0) .and. (dmin <= tol)
                if (ontarget) ontarget = allowed(kbest)
                if ((.not. lists_ok) .or. (.not. ontarget)) then
                    print *, trim(writer) // ': divertor-leg strike ' // &
                        'at ', xs, ys, ' -> nearest wall structure ', &
                        kbest, ' at ', dmin, ' m (on a declared ' // &
                        'target: ', ontarget, ')'
                    if (.not. ontarget) nviol = nviol + 1
                end if
            end do
        end do

        if (nstrike == 0) then
            print *, trim(writer) // ': WARNING: no divertor-leg ' // &
                'strike points found for the given X-point(s)'
        end if

        ! Record problems (do NOT fail here - the caller aggregates and
        ! raises a single error, so the user sees every offending strike)
        nproblems = nviol
        if (.not. lists_ok) then
            print *, trim(writer) // ': not all required divertor ' // &
                'target lists (goat.vessel.solps_*_divertor_*_target) ' // &
                'are given (non-empty)'
            nproblems = max(nproblems, 1_I8)
        elseif (nviol > 0) then
            print *, trim(writer) // ': ', nviol, ' primary strike(s) ' // &
                'are not on any declared divertor target'
        end if

    contains

        ! Number of entries in target list `code`
        function ListCount(code) result(n)
            integer(I8), intent(in) :: code
            integer(I8)             :: n
            select case (code)
            case (1)
                n = int(vessel%nsolpsldi, kind=I8)
            case (2)
                n = int(vessel%nsolpsldo, kind=I8)
            case (3)
                n = int(vessel%nsolpsudi, kind=I8)
            case (4)
                n = int(vessel%nsolpsudo, kind=I8)
            case default
                n = 0
            end select
        end function

        ! Mark the structures of target list `code` in `allowed`
        subroutine MarkList(code)
            integer(I8), intent(in) :: code
            integer(I8)             :: ii
            select case (code)
            case (1)
                do ii = 1, int(vessel%nsolpsldi, kind=I8)
                    allowed(vessel%solpsldiind(ii)) = .true.
                end do
            case (2)
                do ii = 1, int(vessel%nsolpsldo, kind=I8)
                    allowed(vessel%solpsldoind(ii)) = .true.
                end do
            case (3)
                do ii = 1, int(vessel%nsolpsudi, kind=I8)
                    allowed(vessel%solpsudiind(ii)) = .true.
                end do
            case (4)
                do ii = 1, int(vessel%nsolpsudo, kind=I8)
                    allowed(vessel%solpsudoind(ii)) = .true.
                end do
            end select
        end subroutine

    end subroutine

    ! Secondary (additional) X-point inventory and per-X-point divertor
    ! labels. Identifies which X-points are the per-side-innermost
    ! primaries and which are additional (snowflake/secondary nulls),
    ! assigns the per-X divertor label array simgrid%data%xpDivLabel
    ! (0 = primary; else the divertor index the additional X-point's
    ! strike legs hit: SN 1..2 = inner/outer with the USN reversal,
    ! CDN/DDN 1..4 clockwise from the lower inner divertor) and prints a
    ! diagnostic inventory.
    !
    ! NOTE: for now this ONLY does the inventory + labels; it does NOT
    ! yet create the additional volume regions (the wedge extraction is a
    ! later step of the SECONDARY_XPOINT_REGIONS_PLAN). It never raises a
    ! hard error: on any inconsistency it prints a warning and leaves the
    ! label 0.
    subroutine AddSecondaryXPointRegions(simgrid, topomesh, vessel)

        ! Arguments
        type(GridUDT), intent(inout)    :: simgrid
        type(TopomeshUDT), intent(in)   :: topomesh
        type(VesselUDT), intent(in)     :: vessel

        ! Parameters
        real(R8), parameter             :: dtolwall = 5.0e-3_R8

        ! Auxiliary
        integer(I8), allocatable        :: xp(:), sp(:), spx(:), &
            reach(:), wfaces(:), iscompTM(:), wedgeof(:), &
            wcells(:), wptr(:, :), addxi(:), addlab(:), addn(:), &
            newreg(:), order(:), tfc(:), vf(:), bfc(:)
        logical, allocatable            :: isprim(:), istLI(:), &
            istLO(:), istUI(:), istUO(:), istany(:), barrier(:), &
            seed(:), thisw(:), iscore(:), wacc(:)
        logical                         :: okwall, badwall
        real(R8), allocatable           :: addrho(:)
        integer(I8)                     :: nxp, nc, nf, i, j, k, ia, &
            ka, nadd, lab, kbest, nsvi, nregbase, ivadd, wtot, &
            nfixed, nftpfr
        real(R8)                        :: yc, fcore, xs, ys, &
            xxa, yya, xcen, ycen, ux2, uy2, ax, ay, bx, by, &
            wxm, wym, ang, angmin, angmax, cab, caw, cwb
        integer(I8)                     :: sA, sB
        logical                         :: isSN, isupper, ok, okref, okc
        character(64)                   :: lbl, orient

        lbl = simgrid%data%SOLPStopologylabel
        orient = simgrid%data%SOLPStopologyorient
        isSN = (trim(lbl) == 'GEOMETRY_SN')
        isupper = (trim(orient) == 'upper')
        nc = simgrid%cell%ntot
        nf = simgrid%face%ntot

        ! Ensure a valid (zero) label array of the right length.
        ! Narrow-grid rule: X-points without a flux surface in the
        ! delivered grid (orphaned by the narrow-grid target filter,
        ! fieldlineID 0) are dropped BEFORE sizing xpDivLabel, so the
        ! labels stay aligned with the (equally pruned) X-point rows
        ! of the traduit topology block.
        xp = topomesh%GetXPointIDs()
        if (size(xp) > 0) then
            if (any(simgrid%vert%fieldlineID(xp) <= 0)) then
                print *, 'AddSecondaryXPointRegions: ignoring ', &
                    count(simgrid%vert%fieldlineID(xp) <= 0), &
                    ' X-point(s) without a flux surface in the ' // &
                    'filtered grid'
                xp = pack(xp, simgrid%vert%fieldlineID(xp) > 0)
            end if
        end if
        nxp = size(xp)
        if (allocated(simgrid%data%xpDivLabel)) &
            deallocate(simgrid%data%xpDivLabel)
        allocate(simgrid%data%xpDivLabel(nxp))
        simgrid%data%xpDivLabel = 0
        if (nxp == 0) return

        ! Reference + per-side-innermost (primary) X-points
        call SOLPSReferenceCenterFlux(topomesh, yc, fcore, okref)
        if (.not. okref) then
            print *, 'AddSecondaryXPointRegions: no core reference - ' // &
                'all X-points kept primary (label 0)'
            return
        end if
        call ClassifyPrimaryXPoints(topomesh, xp, isprim, okref)
        nadd = nxp - count(isprim)

        print *, 'AddSecondaryXPointRegions: ', nxp, ' X-point(s): ', &
            count(isprim), ' primary (per-side innermost), ', nadd, &
            ' additional'
        if (nadd == 0) then
            print *, 'AddSecondaryXPointRegions: per-X divertor ' // &
                'labels: ', simgrid%data%xpDivLabel
            return
        end if

        ! Target structure masks (divertor-index mapping + wedge wall
        ! validation); istany = the union of all four target lists
        allocate(istLI(vessel%nstructures), istLO(vessel%nstructures), &
            istUI(vessel%nstructures), istUO(vessel%nstructures), &
            istany(vessel%nstructures))
        istLI = .false.
        istLO = .false.
        istUI = .false.
        istUO = .false.
        do i = 1, int(vessel%nsolpsldi, kind=I8)
            istLI(vessel%solpsldiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsldo, kind=I8)
            istLO(vessel%solpsldoind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudi, kind=I8)
            istUI(vessel%solpsudiind(i)) = .true.
        end do
        do i = 1, int(vessel%nsolpsudo, kind=I8)
            istUO(vessel%solpsudoind(i)) = .true.
        end do
        istany = istLI .or. istLO .or. istUI .or. istUO

        sp = topomesh%GetStrikePointIDs()
        spx = topomesh%GetStrikePointXPointIDs()

        ! Core cell mask (from the base cvReg the family writer left:
        ! core = the SOLPS core region IDs 1 (SN) / 1 or 5 (DN)). Used
        ! to sanity-check that a wedge never swallows core cells.
        allocate(iscore(nc))
        do i = 1, nc
            iscore(i) = (simgrid%cell%reg(i) == 1) .or. &
                ((.not. isSN) .and. (simgrid%cell%reg(i) == 5))
        end do
        nregbase = maxval(simgrid%cell%reg)

        ! Per-additional-X arrays and the wedge bookkeeping
        allocate(wedgeof(nc), barrier(nf), seed(nc), thisw(nc), wacc(nc), &
            iscompTM(topomesh%face%ntot), wcells(nc), wptr(nadd, 2), &
            addxi(nadd), addlab(nadd), addn(nadd), addrho(nadd), &
            newreg(nadd), order(nadd))
        wedgeof = 0
        wtot = 0
        ka = 0

        ! ---- Phase 1: per additional X-point: divertor label + wedge ----
        do ia = 1, nxp
            if (isprim(ia)) cycle
            ka = ka + 1
            addxi(ka) = ia
            addrho(ka) = abs(topomesh%vert%fval(xp(ia)) - fcore)

            ! Divertor label = index of the target its strike legs hit
            lab = 0
            nsvi = 0
            do j = 1, size(sp)
                if (spx(j) /= xp(ia)) cycle
                nsvi = nsvi + 1
                xs = topomesh%vert%x(sp(j))
                ys = topomesh%vert%y(sp(j))
                kbest = NearestStructure(xs, ys)
                if (kbest == 0) cycle
                if (lab == 0) lab = DivertorIndex(kbest)
            end do
            simgrid%data%xpDivLabel(ia) = lab
            addlab(ka) = lab

            ! Wedge extraction: barriers = grid faces on this X's
            ! separatrix component; seed the cells at the X grid vertex
            ! on the target side; flood within the barriers.
            call WalkSeparatrixComponent(topomesh, xp(ia), reach, &
                wfaces, ok)
            wptr(ka, 1) = wtot + 1
            wptr(ka, 2) = 0
            addn(ka) = 0
            if (.not. ok) then
                print *, 'AddSecondaryXPointRegions: WARNING: ' // &
                    'separatrix walk failed for additional X-point ', &
                    ia, ' - no wedge region created'
                cycle
            end if
            iscompTM = 0
            iscompTM(wfaces) = 1
            do i = 1, nf
                k = simgrid%face%TMfacelabel(i)
                barrier(i) = .false.
                if ((k >= 1) .and. (k <= topomesh%face%ntot)) then
                    if (iscompTM(k) == 1) barrier(i) = .true.
                end if
            end do

            xxa = topomesh%vert%x(xp(ia))
            yya = topomesh%vert%y(xp(ia))

            ! Wedge extraction (GENERAL -- no dependence on the specific case).
            ! The new PFR of this X is the region its OWN separatrix component
            ! cuts off near the wall. Rather than guessing which strikes are the
            ! real legs (a core-directed arm can strike on the divertor side,
            ! e.g. an SF- X sitting next to a DDN primary X), flood outward from
            ! EVERY grid cell touching the X grid vertex, blocked by the
            ! component. Each sector at the X floods either into a small ENCLOSED
            ! lobe piece (kept) or into the big SOL/core rest (which reaches core
            ! cells or most of the grid, and is rejected). The union of the kept
            ! pieces is the wedge; this also merges the two halves of a wedge that
            ! the central continuation over-splits. Works for any topology and any
            ! number/position of the X's strikes.
            ivadd = NearestGridVertex(xxa, yya)
            vf = GetVertFace(simgrid%vert, ivadd)
            wacc = .false.
            do j = 1, size(vf)
                tfc = GetFaceCell(simgrid%face, vf(j))
                do k = 1, size(tfc)
                    if (wacc(tfc(k))) cycle
                    seed = .false.
                    seed(tfc(k)) = .true.
                    call FloodWedge()
                    if (count(thisw) == 0) cycle
                    if (any(thisw .and. iscore)) cycle    ! big SOL/core rest
                    if (count(thisw) > nc/2) cycle         ! runaway (rest)
                    ! A genuine private lobe touches the wall ONLY on declared
                    ! target structures (its mouth is target wall by the
                    ! preprocessing contract). A SOL pocket sealed by the
                    ! separatrix arms can also be non-core and < nc/2, but it
                    ! touches the MAIN (non-target) wall -> reject it, so only the
                    ! true lobe(s) survive. Require at least one (target) wall
                    ! face too, so a fully-interior pocket is not taken.
                    okwall = .false.
                    badwall = .false.
                    do i = 1, nf
                        if (.not. simgrid%face%BF(i)) cycle
                        bfc = GetFaceCell(simgrid%face, i)
                        if (.not. thisw(bfc(1))) cycle
                        okwall = .true.
                        xs = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) + &
                                     simgrid%vert%x(simgrid%face%vert(i, 2)))
                        ys = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) + &
                                     simgrid%vert%y(simgrid%face%vert(i, 2)))
                        kbest = NearestStructure(xs, ys)
                        if ((kbest == 0) .or. (.not. istany(kbest))) then
                            badwall = .true.
                            exit
                        end if
                    end do
                    if ((.not. okwall) .or. badwall) cycle
                    wacc = wacc .or. thisw                 ! an enclosed lobe piece
                end do
            end do
            thisw = wacc

            if (count(thisw) == 0) then
                print *, 'AddSecondaryXPointRegions: WARNING: no ' // &
                    'enclosed wedge found for additional X-point ', ia, &
                    ' - no wedge region created'
                cycle
            end if
            do i = 1, nc
                if (thisw(i)) then
                    wtot = wtot + 1
                    wcells(wtot) = i
                end if
            end do
            wptr(ka, 2) = wtot - wptr(ka, 1) + 1
            addn(ka) = wptr(ka, 2)
            print *, 'AddSecondaryXPointRegions: additional X-point ', &
                ia, ' at ', xxa, yya, ' label ', lab, ' -> wedge of ', &
                addn(ka), ' cells'
        end do

        ! ---- Phase 2: assign wedgeof (largest first, so a nested ----
        ! inner wedge, assigned last, wins its shared cells) ----
        do i = 1, nadd
            order(i) = i
        end do
        call SortWedges()          ! order() by addn descending
        do i = 1, nadd
            ka = order(i)
            do j = wptr(ka, 1), wptr(ka, 1) + wptr(ka, 2) - 1
                wedgeof(wcells(j)) = ka
            end do
        end do

        ! ---- Phase 3: region numbering (spec 2.9): additional regions ----
        ! appended after the basic NREG, ordered by (label asc, rho asc,
        ! outer-first). Compute a rank per wedge. ----
        do i = 1, nadd
            order(i) = i
        end do
        call SortNumbering()       ! order() by (label, rho, -n)
        newreg = 0
        k = 0
        do i = 1, nadd
            ka = order(i)
            if (addn(ka) == 0) cycle       ! no valid wedge -> no region
            k = k + 1
            newreg(ka) = nregbase + k
        end do

        ! ---- Phase 4: overwrite cvReg of the wedge cells ----
        do i = 1, nc
            if (wedgeof(i) > 0) then
                if (newreg(wedgeof(i)) > 0) &
                    simgrid%cell%reg(i) = newreg(wedgeof(i))
            end if
        end do

        ! ---- Phase 5: face-region patches for wedge cells ----
        ! (a) every boundary face adjacent to a wedge cell must be on a
        !     declared target (spec 3.2 fail-fast);
        ! (b) chain faces that end up interior to a wedge -> fcReg 0.
        nfixed = 0
        do i = 1, nf
            if (simgrid%face%BF(i)) then
                tfc = GetFaceCell(simgrid%face, i)
                if (wedgeof(tfc(1)) == 0) cycle
                ! wedge-adjacent boundary face: require a target wall
                xs = 0.5_R8*(simgrid%vert%x(simgrid%face%vert(i, 1)) + &
                             simgrid%vert%x(simgrid%face%vert(i, 2)))
                ys = 0.5_R8*(simgrid%vert%y(simgrid%face%vert(i, 1)) + &
                             simgrid%vert%y(simgrid%face%vert(i, 2)))
                kbest = NearestStructure(xs, ys)
                if (kbest == 0) then
                    call gdErrorHandler('AddSecondaryXPointRegions: a ' // &
                        'wedge (additional-PFR) boundary face is not ' // &
                        'on any wall structure - extend the ' // &
                        'goat.vessel.solps_*_target lists so the whole ' // &
                        'wall span of every additional X-point is ' // &
                        'declared target wall')
                elseif (.not. istany(kbest)) then
                    call gdErrorHandler('AddSecondaryXPointRegions: a ' // &
                        'wedge (additional-PFR) boundary face lands on ' // &
                        'a non-target wall structure - extend the ' // &
                        'goat.vessel.solps_*_target lists to cover it')
                end if
            else
                ! interior face fully inside a wedge that carried a chain
                ! (entrance/cut) region -> reset to 0 ("split in two")
                tfc = GetFaceCell(simgrid%face, i)
                if (size(tfc) /= 2) cycle
                if ((wedgeof(tfc(1)) > 0) .and. (wedgeof(tfc(2)) > 0) &
                    .and. (simgrid%face%reg(i) /= 0)) then
                    simgrid%face%reg(i) = 0
                    nfixed = nfixed + 1
                end if
            end if
        end do
        if (nfixed > 0) print *, 'AddSecondaryXPointRegions: reset ', &
            nfixed, ' interior chain faces inside wedges to fcReg 0'

        ! ---- Phase 6: flux tubes fully inside a wedge -> ftReg 3 (PFR) ----
        nftpfr = 0
        associate(fd => simgrid%data%fluxdata)
        do i = 1, fd%nFt
            ok = .true.
            do j = fd%fluxtubecellsP(i, 1), &
                fd%fluxtubecellsP(i, 1) + fd%fluxtubecellsP(i, 2) - 1
                if (wedgeof(fd%fluxtubecells(j)) == 0) then
                    ok = .false.
                    exit
                end if
            end do
            if (ok) then
                fd%fluxtuberegID(i) = 3
                nftpfr = nftpfr + 1
            end if
        end do
        end associate

        print *, 'AddSecondaryXPointRegions: created ', k, &
            ' additional volume region(s) ', nregbase + 1, '..', &
            nregbase + k, '; ', nftpfr, ' flux tubes set to PFR'
        print *, 'AddSecondaryXPointRegions: per-X divertor labels: ', &
            simgrid%data%xpDivLabel

    contains

        ! Divertor index of structure kk: CDN/DDN 1..4 (LI/UI/UO/LO,
        ! clockwise from the lower inner target), SN 1..2 (inner/outer,
        ! reversed for an upper single null)
        function DivertorIndex(kk) result(idx)
            integer(I8), intent(in) :: kk
            integer(I8)             :: idx
            idx = 0
            if (isSN) then
                if (isupper) then
                    if (istUO(kk)) idx = 1
                    if (istUI(kk)) idx = 2
                else
                    if (istLI(kk)) idx = 1
                    if (istLO(kk)) idx = 2
                end if
            else
                if (istLI(kk)) idx = 1
                if (istUI(kk)) idx = 2
                if (istUO(kk)) idx = 3
                if (istLO(kk)) idx = 4
            end if
        end function

        ! Nearest wall sub-structure to (px,py) within dtolwall, else 0
        function NearestStructure(px, py) result(ks)
            real(R8), intent(in)    :: px, py
            integer(I8)             :: ks, kk
            real(R8)                :: dd, dm
            ks = 0
            dm = posinfval_R8()
            do kk = 1, int(vessel%nstructures, kind=I8)
                dd = StructureDistance(vessel, px, py, kk)
                if (dd < dm) then
                    dm = dd
                    ks = kk
                end if
            end do
            if (dm > dtolwall) ks = 0
        end function

        ! Nearest grid vertex to (px,py)
        function NearestGridVertex(px, py) result(iv)
            real(R8), intent(in)    :: px, py
            integer(I8)             :: iv, ii
            real(R8)                :: dd, dm
            iv = 1
            dm = posinfval_R8()
            do ii = 1, simgrid%vert%ntot
                dd = (simgrid%vert%x(ii) - px)**2 + &
                     (simgrid%vert%y(ii) - py)**2
                if (dd < dm) then
                    dm = dd
                    iv = ii
                end if
            end do
        end function

        ! Flood the wedge into thisw() from seed(), blocked by barrier()
        ! faces (host-associated); wall (boundary) faces stop the fill.
        subroutine FloodWedge()
            integer(I8), allocatable    :: q(:), cfl(:), tfcl(:)
            integer(I8)                 :: h, t, cq, jj, kk, cn
            allocate(q(nc))
            thisw = .false.
            h = 0
            t = 0
            do jj = 1, nc
                if (seed(jj) .and. (.not. thisw(jj))) then
                    t = t + 1
                    q(t) = jj
                    thisw(jj) = .true.
                end if
            end do
            h = 0
            do while (h < t)
                if (t > nc/2) exit     ! runaway: the big SOL/core rest, not a
                                       ! lobe -- the caller rejects it on size
                h = h + 1
                cq = q(h)
                cfl = GetCellFace(simgrid%cell, cq)
                do jj = 1, size(cfl)
                    if (simgrid%face%BF(cfl(jj))) cycle
                    if (barrier(cfl(jj))) cycle
                    tfcl = GetFaceCell(simgrid%face, cfl(jj))
                    do kk = 1, size(tfcl)
                        cn = tfcl(kk)
                        if (cn == cq) cycle
                        if (thisw(cn)) cycle
                        thisw(cn) = .true.
                        t = t + 1
                        q(t) = cn
                    end do
                end do
            end do
        end subroutine

        ! order() <- wedge indices sorted by cell count descending
        subroutine SortWedges()
            integer(I8) :: a, b, tmp
            do a = 1, nadd - 1
                do b = a + 1, nadd
                    if (addn(order(b)) > addn(order(a))) then
                        tmp = order(a); order(a) = order(b); order(b) = tmp
                    end if
                end do
            end do
        end subroutine

        ! order() <- wedge indices sorted by (label asc, rho asc, n desc)
        subroutine SortNumbering()
            integer(I8) :: a, b, tmp
            logical     :: swap
            do a = 1, nadd - 1
                do b = a + 1, nadd
                    swap = .false.
                    if (addlab(order(b)) < addlab(order(a))) then
                        swap = .true.
                    elseif (addlab(order(b)) == addlab(order(a))) then
                        if (addrho(order(b)) < addrho(order(a)) - 1.0e-12_R8) then
                            swap = .true.
                        elseif (abs(addrho(order(b)) - addrho(order(a))) &
                                <= 1.0e-12_R8) then
                            if (addn(order(b)) > addn(order(a))) swap = .true.
                        end if
                    end if
                    if (swap) then
                        tmp = order(a); order(a) = order(b); order(b) = tmp
                    end if
                end do
            end do
        end subroutine

    end subroutine

    ! Reference center height and core-side flux value (as in
    ! ClassifyBasicSOLPSCatalogTopology): from the O-point vertex when
    ! the axis is meshed, otherwise from the core-cell centroid and the
    ! core boundary faces. ok=.false. if there is no usable reference.
    subroutine SOLPSReferenceCenterFlux(topomesh, yc, fcore, ok)
        type(TopomeshUDT), intent(in)   :: topomesh
        real(R8), intent(out)           :: yc, fcore
        logical, intent(out)            :: ok
        integer(I8), allocatable        :: op(:), cc(:), ccf(:), tmpi(:)
        integer(I8)                     :: i
        ok = .true.
        yc = 0.0_R8
        fcore = 0.0_R8
        op = topomesh%GetOPointIDs()
        if (size(op) == 1) then
            yc = topomesh%vert%y(op(1))
            fcore = topomesh%vert%fval(op(1))
            return
        end if
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            ok = .false.
            return
        end if
        allocate(tmpi(0))
        do i = 1, size(cc)
            tmpi = [tmpi, topomesh%cell%GetVert(cc(i))]
        end do
        yc = sum(topomesh%vert%y(tmpi))/real(size(tmpi), kind=R8)
        deallocate(tmpi)
        ccf = topomesh%GetCoreFaceIDs()
        if (size(ccf) == 0) then
            ok = .false.
            return
        end if
        do i = 1, size(ccf)
            fcore = fcore + 0.5_R8*( &
                topomesh%vert%fval(topomesh%face%vert(ccf(i), 1)) + &
                topomesh%vert%fval(topomesh%face%vert(ccf(i), 2)))
        end do
        fcore = fcore/real(size(ccf), kind=R8)
    end subroutine

    ! Core center point (R,Z): the O-point vertex when the axis is
    ! meshed, otherwise the core-cell vertex centroid. ok=.false. when
    ! there is no core. Used to tell a divertor leg (pointing away from
    ! the core) from a secondary-separatrix grazing arm.
    subroutine SOLPSCoreCenter(topomesh, xc, yc, ok)
        type(TopomeshUDT), intent(in)   :: topomesh
        real(R8), intent(out)           :: xc, yc
        logical, intent(out)            :: ok
        integer(I8), allocatable        :: op(:), cc(:), tmpi(:)
        integer(I8)                     :: i
        ok = .true.
        xc = 0.0_R8
        yc = 0.0_R8
        op = topomesh%GetOPointIDs()
        if (size(op) == 1) then
            xc = topomesh%vert%x(op(1))
            yc = topomesh%vert%y(op(1))
            return
        end if
        cc = topomesh%GetCoreCellIDs()
        if (size(cc) == 0) then
            ok = .false.
            return
        end if
        allocate(tmpi(0))
        do i = 1, size(cc)
            tmpi = [tmpi, topomesh%cell%GetVert(cc(i))]
        end do
        xc = sum(topomesh%vert%x(tmpi))/real(size(tmpi), kind=R8)
        yc = sum(topomesh%vert%y(tmpi))/real(size(tmpi), kind=R8)
    end subroutine

    ! Flag the per-side-innermost (primary) X-points among the topomesh
    ! saddles xp = topomesh%GetXPointIDs(): the one nearest in flux to
    ! the core reference on each side (below/above yc). isprim is sized
    ! like xp. ok=.false. if there is no usable reference.
    subroutine ClassifyPrimaryXPoints(topomesh, xp, isprim, ok)
        type(TopomeshUDT), intent(in)           :: topomesh
        integer(I8), intent(in)                 :: xp(:)
        logical, allocatable, intent(out)       :: isprim(:)
        logical, intent(out)                    :: ok
        real(R8)                                :: yc, fcore, d, dlow, dup
        integer(I8)                             :: i, ilow, iup
        allocate(isprim(size(xp)))
        isprim = .false.
        call SOLPSReferenceCenterFlux(topomesh, yc, fcore, ok)
        if (.not. ok) return
        ilow = 0
        iup = 0
        dlow = 0.0_R8
        dup = 0.0_R8
        do i = 1, size(xp)
            d = abs(topomesh%vert%fval(xp(i)) - fcore)
            if (topomesh%vert%y(xp(i)) >= yc) then
                if ((iup == 0) .or. (d < dup)) then
                    iup = i
                    dup = d
                end if
            else
                if ((ilow == 0) .or. (d < dlow)) then
                    ilow = i
                    dlow = d
                end if
            end if
        end do
        if (ilow > 0) isprim(ilow) = .true.
        if (iup > 0) isprim(iup) = .true.
    end subroutine

    ! Minimum distance from a point to vessel sub-structure ks
    function StructureDistance(vessel, px, py, ks) result(dres)
        type(VesselUDT), intent(in) :: vessel
        real(R8), intent(in)        :: px, py
        integer(I8), intent(in)     :: ks
        real(R8)                    :: dres
        integer(I8)                 :: jj, npt
        real(R8)                    :: x1, y1, x2, y2
        dres = posinfval_R8()
        npt = vessel%structures(ks)%np
        do jj = 1, npt - 1
            x1 = vessel%structures(ks)%x(jj)
            y1 = vessel%structures(ks)%y(jj)
            x2 = vessel%structures(ks)%x(jj + 1)
            y2 = vessel%structures(ks)%y(jj + 1)
            dres = min(dres, SegmentDistance(px, py, x1, y1, x2, y2))
        end do
        if (vessel%structures(ks)%isclosed .and. (npt > 2)) then
            x1 = vessel%structures(ks)%x(npt)
            y1 = vessel%structures(ks)%y(npt)
            x2 = vessel%structures(ks)%x(1)
            y2 = vessel%structures(ks)%y(1)
            dres = min(dres, SegmentDistance(px, py, x1, y1, x2, y2))
        end if
    end function

    ! Minimum distance from a point to one segment
    function SegmentDistance(px, py, x1, y1, x2, y2) result(dres)
        real(R8), intent(in)    :: px, py, x1, y1, x2, y2
        real(R8)                :: dres
        real(R8)                :: dx, dy, t
        dx = x2 - x1
        dy = y2 - y1
        if ((dx*dx + dy*dy) <= 0.0_R8) then
            dres = sqrt((px - x1)**2 + (py - y1)**2)
            return
        end if
        t = ((px - x1)*dx + (py - y1)*dy)/(dx*dx + dy*dy)
        t = max(0.0_R8, min(1.0_R8, t))
        dres = sqrt((px - (x1 + t*dx))**2 + (py - (y1 + t*dy))**2)
    end function


end module
