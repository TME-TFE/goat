!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General definitions to be used in other modules. These definitions can
! be for example the integer that is linked to a specific type of 
! vessel boundary (e.g. target plate). These definitions should replace
! hard-coded constants etc. 

module mod_definitions

    use mod_precision

    implicit none
    public 

    ! Goat
    !=====
    ! Version
    character(*), parameter :: goatversion = '01.000.000s'

    ! Definitions for grid deformation
    !=================================
    ! Vessel parts definitions
    integer(I8), parameter  :: targetID = 1, coreID = 2, outerboundaryID = 3, &
        vesselID = 4, interiorID = 5
    
    ! Definitions for grid generation
    !================================
    ! Topological mesh vertex IDs (field minimum, saddle point, maximum: 1, 2, 3,
    ! field & boundary tangency point types 1 and 2: 4, 5, boundary vertex
    ! that is no tangency point or whatever: 6, vertex constructed by 
    ! splitting face: -1)
    ! Note: type 1 TP contour line goes out of domain, type 2 goes inwards
    integer(I8), parameter  :: TMvertexminID = 1, TMvertexsaddleID = 2, &
        TMvertexmaxID = 3, TMvertextp1ID = 4, TMvertextp2ID = 5, TMvertexbndID = 6, &
        TMvertexsplitID = -1

    ! Topological mesh boundary IDs (1: radial face, 2: poloidal face, &
    ! 3: non-aligned boundary face, 4: separatrix/saddle point, 5: inserted core boundary, &
    ! 6: inserted PF boundary), 7: aligned boundary face, 8: scrape-off layer 
    ! flux surface (aligned boundary) ID
    integer(I8), parameter :: TMfaceradID = 1, TMfacepolID = 2, TMfacebndID = 3, &
        TMfacesepID = 4, TMfacecoreID = 5, TMfacePFID = 6, TMfacealbndID = 7, &
        TMfaceSOLID = 8

    ! Joint IDs: aligned faces (TMfacealignedID), non-aligned faces, 
    ! domain boundary faces (a bit of a sad choice of names has led to this)
    integer(I8), parameter, dimension(*) :: &
        TMfacealignedID = [TMfacepolID, TMfacesepID, TMfacecoreID, TMfacePFID, TMfacealbndID, TMfaceSOLID], &
        TMfacenonalignedID = [TMfaceradID, TMfacebndID], &
        TMfaceBFID = [TMfacebndID, TMfacealbndID, TMfacecoreID, TMfacePFID, TMfaceSOLID]

    ! Definitions for SOLPS
    !======================
    ! Version
    character(*), parameter :: SOLPSversion = '03.002.000'

    ! Solver-related
    integer(I8), parameter :: SOLPScoreregID = 1, SOLPScoreregIDincr = 4, &
        SOLPSbndcellID = 3, SOLPSinternalcellID = 1, SOLPScorefclblID = -21, &
        SOLPScorefclblIDincr = -4
   
    ! Topological mesh identification number
    integer(I8), parameter :: TMTopSN = 1, TMTopDN = 2, TMTopGeneral = 0, &
        TMTopL = 3, TMTopLM = 4
end module