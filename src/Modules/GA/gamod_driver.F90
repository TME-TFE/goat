!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all routines related to the driver for the grid adaptation module, such as initialization and postprocessing routines.

module gamod_driver

    ! Initialize
    !============
    ! Load modules
    use mod_triangulation
    use UnstructuredInterpolant2D
    use goatmod_types
    use goatmod_userinput
    use gamod_utility
    use gamod_types

    ! The usual
    implicit none
    save
    public     
    
    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!  
    
    contains

    subroutine GridAdaptor(grid, environment, magneticField, options)

        ! Description
        !============
        ! Overarching driver for grid adaptation, one lever lower than GADriver. 
        ! Adapts topology of unstructured grid with the goal to improve grid quality based on grid metric and user inputs.

        ! Initialize
        !===========
        ! Declare modules  
        
        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)              :: grid
        type(EnvironmentUDT), intent(in)            :: environment
        type(MagneticFieldUDT), intent(in)          :: magneticField     
        type(GAoptionsUDT), intent(inout)           :: options

        ! Initialize grid adaptation
        !===========================
        call GAinit(grid, options, environment, magneticField)
        

        ! Driver Selection
        !=================
        select case (options%meth)

        case ('simple')

            ! Regular grid adaption
            call GAInternalDriver(grid, options, environment, magneticField)

        case ('aposteriori')

            ! Grid adaptation based on simulation information
            call GAapostDriver(grid, options, environment, magneticField) 

        case default

            ! Call error handler
            call gdErrorHandler('GridAdaptor: unknown driver option')
        
        end select

        ! Postprocessing
        !===============
        call PostProcessGA(grid, magneticField, options)

    end subroutine

    subroutine GAinit(grid,options,environment,magneticField)

        ! Description
        !============

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)          :: grid 
        type(GAoptionsUDT), intent(inout)       :: options
        type(EnvironmentUDT), intent(in)        :: environment
        type(MagneticFieldUDT), intent(in)      :: magneticField

        ! Variables
        integer(I8) :: i
        integer(I8), allocatable, dimension(:) :: tv, cvLookUp, indcv
        logical :: cells(grid%cell%ntot), is_ordered(grid%cell%ntot), &
            use_nsep, use_sepID, start, b_flag(grid%face%ntot)
        character(:), allocatable :: base_func


        ! Initialize
        !===========
 
        associate(&
            c  => grid%cell, &
            v  => grid%vert &
            )

        ! Recompute cell centers
        do i = 1, c%ntot
            ! Get cell vertices
            tv = GetCellVertGA(c, i)

            ! Compute coordinates
            call c%x%SetSingleElement(i,sum(v%x%Get(tv))/real(size(tv), kind=R8)) 
            call c%y%SetSingleElement(i,sum(v%y%Get(tv))/real(size(tv), kind=R8)) 

        end do

        ! Determine cflags
        b_flag = [grid%face%label%Get() /= 0]
        indcv = (/(i, i = 1, c%ntot)/)
        call grid%DetermineCflags(indcv, b_flag)


        ! Check order of vertices  (see GetGeo_usCouples.m)
        cells = .true.
        call grid%CheckVertOrder(is_ordered, cells) 

        if (.not. all(is_ordered)) then

            ! Do ordening
            call grid%ReorderCellConn(is_ordered)

        end if 

        ! Check fsFc
        call grid%CheckFsFc()

        ! Get fsVx from fsFc
        call grid%GetFsVxFromFsFc(options)

        ! Determine Xpoints and separatrices
        cvLookUp = GetCvLookUpGA(c)
        use_nsep = .false.
        call grid%GiveXpoints(use_nsep,cvLookUp)
        use_sepID = .false.
        start = .true.
        call grid%GiveSeparatrices(use_nsep, use_sepID, start, cvLookUp)

        ! Identify aligned faces
        call grid%IdentifyAlignedFaces(options,magneticField)

        ! Set up the distance functions
        if (options%dist_function) then

            ! Define base function to compute the distance function
            base_func = options%base_func !'exp(-dist/d)'

            ! Free distance function based on user-defined input
            grid%fun%d_char_type = options%d_char_type
            grid%fun%dist_type   = options%dist_type
            grid%fun%d_rescale   = options%d_rescale
            call grid%fun%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for high poloidal flux next to the separatrix
            grid%fun_r%d_char_type = options%d_char_type
            grid%fun_r%dist_type   = 'pol_flux_est'
            grid%fun_r%d_rescale   = options%d_rescale
            call grid%fun_r%ComputeDistanceFunction(grid,options,base_func)

            ! Distance function for wall proximity
            grid%fun_wall%d_char_type = options%d_char_type
            grid%fun_wall%dist_type = options%dist_type_wall
            grid%fun_wall%d_rescale = options%d_rescale_wall
            call grid%fun_wall%ComputeDistanceFunction(grid,options,base_func)

        end if 


        ! Check the connectivity
        if (options%debug) then
            call grid%CheckUnstructuredGrid()
        end if

        ! Correct face labels on for wide grid
        call grid%CheckFcLbl(options)

        ! Check consistency of options before inspecting the option arrays
        call CheckGAoptions(options)

        ! Construct the far-SOL selector only when an enabled adaptation
        ! criterion actually consumes it. Unconditional construction makes
        ! unrelated operations fail on otherwise valid general topologies.
        if (UsesFarSOLCriteria(options)) then
            call grid%IdentifyfarSOLcells(options)
        end if

        ! Visualize starting grid
        call grid%WriteData('grid_before_GA')
        call grid%WriteFluxSurfaceData()

        end associate        


    end subroutine

    logical function UsesFarSOLCriteria(options) result(usesFarSOL)

        ! Return whether any enabled splitting or merging pass requires the
        ! far-SOL interpolant constructed by IdentifyfarSOLcells.

        type(GAoptionsUDT), intent(in) :: options

        integer(I8) :: i

        usesFarSOL = .false.

        do i = 1, size(options%merging_array)
            if (options%merging_array(i) == 1 .and. &
                options%merge_crit_array(i) == 8) then
                usesFarSOL = .true.
                return
            end if
        end do

        do i = 1, size(options%splitting_array)
            if (options%splitting_array(i) /= 1) cycle
            if (options%splittype_array(i) == 1 .and. &
                any(options%rad_type_array(i) == [5, 6])) then
                usesFarSOL = .true.
                return
            end if
            if (options%splittype_array(i) == 2 .and. &
                any(options%pol_type_array(i) == [4, 5])) then
                usesFarSOL = .true.
                return
            end if
        end do

    end function

    subroutine GAInternalDriver(grid,options,environment,magneticField)

        ! Description
        !============
        ! Internal driver for the grid adaptation where all real adaptation take place such as removal of small triangles, stacked triangles

        ! Declare variables
        !==================
        ! Argument
        type(GAGridUDT), intent(inout)       :: grid
        type(GAoptionsUDT), intent(inout)    :: options
        type(EnvironmentUDT), intent(in)     :: environment
        type(MagneticFieldUDT), intent(in)   :: magneticField

        ! Auxiliary
        integer(I8) :: i
        type(QualityMetricUDT) :: qm
        type(GAoptionsUDT) :: options1
        type(GAoptionsUDT) :: options_merge
        type(GAoptionsUDT) :: options_split

        ! Calculate quality metric
        call qm%Initialize(grid)
        call qm%CalculateQualityMetrics(grid, options, magneticField,.false.,.false.)

        ! Remove Small triangles
        if (options%rem_small_trias) &
            call grid%RemoveSmallTriangle(magneticField, qm, options)

        ! Visualize starting grid
        call grid%WriteData('grid_after_rem_trias')

        ! Remove flux tubes with only two triangles
        if (options%rem_trias_tube .or. options%rem_outershell) &
            call grid%RemTriasFlux(options)

        ! Stacked to cutcell
        if (options%stacked_to_cutcell) &
            call grid%StackedToCutcell(magneticField, options)

        ! Splitting non-aligned quads
        if (options%split_noalignedquads) then
            options1 = options
            options1%splittype = 'rad'
            options1%rad_type = 7 !'no_aligned_faces'
            options1%n_split = grid%cell%ntot
            call grid%DoSplitting(magneticField, qm, options1)
        end if

        ! Splitting trapezoids in concave shaved-off flux tube
        if (options%split_shaved_off_tube) then
            options1 = options
            options1%splittype = 'rad'
            options1%rad_type = 8 !'shaved-off_tubes'
            options1%n_split = grid%cell%ntot
            call grid%DoSplitting(magneticField, qm, options1)
        end if

        ! Splitting  and merging
        ! Merging
        do i = 1, size(options%merge_crit_array)
            if (options%merging_array(i) == 1) then
                options_merge = options
                options_merge%merging = .true.
                options_merge%n_merge = options%n_merge_array(i)
                options_merge%merge_crit = options%merge_crit_array(i)
                call grid%DoMerging(magneticField, qm, options_merge)
            end if

            if (options%splitting_array(i) == 1) then
                options_split = options
                options_split%splitting = .true.
                options_split%n_split = options%n_split_array(i)
                options_split%rad_type = options%rad_type_array(i)
                options_split%pol_type = options%pol_type_array(i)
                if (options%splittype_array(i) == 1) then
                    options_split%splittype = 'rad'
                else if (options%splittype_array(i) == 2) then
                    options_split%splittype = 'pol'
                end if
                call grid%DoSplitting(magneticField, qm, options_split)
            end if
        end do

        ! Stacked triangles
        if (options%stacked_trias) &
            call grid%StackedTrias(magneticField, qm, options)

        ! Remove sticking out triangles
        if (options%rem_stickout_trias) &
            call grid%RemoveStickOutTrias(options)

        ! Remove boundary flux tubes with only two triangles
        if (options%rem_trias_tube .or. options%rem_outershell) &
            call grid%RemTriasFlux(options)

        ! Remove stickout quad
        if (options%rem_stickout_quad) &
            call grid%RemoveStickoutQuads()

        ! Boundary layer grid
        if (options%BLG) then
            call grid%BoundaryLayerGrid(qm, magneticField, options)
        end if


    end subroutine

    subroutine GAapostDriver(grid, options, environment, magneticField)

        ! Description
        !============
        ! Internal driver of the grid adaptation for refinement where the computed
        ! residuals or state gradients are high

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)      :: grid
        type(GAoptionsUDT), intent(inout)   :: options
        type(EnvironmentUDT), intent(in)    :: environment
        type(MagneticFieldUDT), intent(in)  :: magneticField

        ! Auxiliary
        integer(I8) :: j, split_cv
        type(TriangulationUDT)              :: triangulation
        type(StateUDT)                      :: state
        type(StateUDT)                      :: state_v
        type(StateUDT)                      :: state_int
        type(UnstructuredInterpolant2DUDT)  :: interp
        type(QualityMetricUDT)              :: qm

        ! Unpack state
        state = environment%SOLPSstate

        ! Interpolation
        !==============
        ! Triangulate
        call grid%TriangulateGAGrid(triangulation)

        ! Interpolate state to vertex positions
        call grid%InterpolateCvToVx(options, state, state_v)

        ! Construct interpolant
        select case (options%apost_interpolation_meth)
        case ('barycentric')

            call interp%SetParametersUS(options%apost_interpolation_meth, 0, 1, triangulation, 'polynomial')

        case ('finite_element')

            call interp%SetParametersUS(options%apost_interpolation_meth, &
                options%apost_interpolationC, options%apost_interpolationM, triangulation, 'polynomial')

        end select



        ! Convert stacked triangle back to cutcells
        if (options%stacked_to_cutcell) &
            call grid%StackedToCutcell(magneticField, options)
        
        ! Splitting
        print *, 'Starting aposteriori splitting'

        ! Select cell to split
        call grid%SelectSplitCellAposteriori(magneticField, options, interp, state, split_cv)
        call qm%Initialize(grid)

        ! While a splitting cell is found
        j = 0
        do while (split_cv /= 0 .and. j .lt. options%n_split)

            if (options%debug) print *, 'Cell: ', split_cv

            call grid%Splitting(split_cv, options, magneticField)

            if (.not.options%slab) call grid%RecalcMagn(magneticField)

            ! Update counter
            j = j + 1

            ! Printing 
            print *, 'Aposteriori splitting: ', j

            ! Check grid
            if (options%debug) call grid%CheckUnstructuredGrid()

            if (options%splittype == 'rad') then

                ! Remove small triangle which were possibly created
                ! Calculate quality metric
                call qm%CalculateQualityMetrics(grid, options, magneticField, .false., .false.)

                ! Remove Small triangles
                if (options%rem_small_trias) then

                    print *, 'Removing small triangles'
                    call grid%RemoveSmallTriangle(magneticField, qm, options)
                    print *, 'Ended removing small triangles'

                    if (options%debug) call grid%CheckUnstructuredGrid()

                end if

            end if

            if (options%apost_meth == 'res') then
                print *, 'Warning: here should be a routine to reduce the residuals of ' // &
                ' the splitting cell to avoid resplitting the same cell over and over.' //&
                ' In this implementation state_v should be changed '// &
                'as this is the state used to interpolate from.'
            end if

            ! Interpolant state on new grid (state = state_int)
            call grid%InterpolateState(interp, state_v, options, state_int)

            ! Select cell to split
            call grid%SelectSplitCellAposteriori(magneticField, options, interp, state_int, split_cv)

            
        end do

        ! Transform remaining pentagons into triangles
        if (options%pents_to_tria) call grid%TransPentsToTrias()

        ! Ordening not necessary probably

        ! Display progress
        print *, 'Ended aposteriori splitting'


    end subroutine
 
    subroutine PostProcessGA(grid,magneticField,options)

        ! Description
        !============
        ! Postprocessing specifically after grid adaptations

        ! Declare variables
        !==================
        ! Arguments
        type(GAGridUDT), intent(inout)      :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(in)         :: options

        ! Auxiliary
        logical :: is_ordered(grid%cell%ntot), cells(grid%cell%ntot), &
         use_nsep, use_sepID, start, err
        integer(I8), allocatable :: cvLookUp(:), fcs(:), f_ord(:,:), nf(:), &
            lbls(:), lbls2(:), fsVx(:), verts(:), fcsv(:), ind(:)
        integer(I8) ::  i, iv, nl, nvi, lb, fcReg(grid%face%ntot), &
            fcLbl_loc(grid%face%ntot), indFc(grid%face%ntot), &
            nflbl
        

        ! Check consistency
        if (options%debug) then
            call grid%CheckUnstructuredGrid()
        end if

        ! Order correctly
        cells = .true.
        call grid%CheckVertOrder(is_ordered,cells) 
        call grid%ReorderCellConn(is_ordered)

        ! Determine number of Xpoints and separatrices
        cvLookUp = GetCvLookUpGA(grid%cell)
        use_nsep = .true.
        use_sepID = .true. 
        call grid%GiveSeparatrices(use_nsep, use_sepID, start, cvLookUp)
        call grid%GiveXpoints(use_sepID,cvLookUp)

        ! Check flux surfaces for possible merging of open surfaces
        call grid%MergeFS()

        ! Recalculate bx, by
        call grid%RecalcMagn(magneticField)

        ! Check whether every vertices is an flux surface for the correct cases
        if (options%stacked_trias .and. .not.options%vesselmode) then

            ! Get all vertices belonging to a flux surface
            fsVx = grid%data%fluxdata%fluxsurfaceverts%Get()

            do iv = 1, grid%vert%ntot

                nvi = count(fsVx == iv)

                if (nvi /= 1) then

                    ! Only allowed when vertex is a boundary vertex and 
                    ! its boundary faces are not aligned
                    err = .false.

                    if (.not.isBoundaryVertGA(grid, iv)) then

                        err = .true.

                    else 

                        ! Get faces of vertices
                        fcsv = GetVertFaceGA(grid%face, iv)
                        if (count(grid%face%aligned%Get(fcsv) == 1) .gt. 0) err = .true.

                    end if

                    if (err) then

                        ! Give error information
                        print *, 'Vertex without flux surface: ', iv
                        print *, grid%vert%x%Get(iv)
                        print *, grid%vert%y%Get(iv)
                        verts = [iv, iv]
                        call grid%WriteErrorData(verts, 1)
                        call gdErrorHandler('PostprocessGA: vertex does not occur once in fsVx')

                    end if

                end if 

            end do

        end if

        ! fcReg
        !------
        ! inner main target => fcReg == 1 
        ! outer main target => fcReg == 4 
        ! inner secondary target => fcReg == 5
        ! inner secondary target => fcReg == 8  
        
        ! Locally 
        fcLbL_loc = GetfcLblGA(grid%face,options)

        call Unique(fcLbl_loc, lbls)
        allocate(lbls2(count(lbls /= 0)))
        lbls2 = pack(lbls,lbls /= 0)
        nl = size(lbls2)

        if (nl .gt. size(options%fcRegmappingGA)) then
            call gdErrorHandler('PostProcesGA: fcReg mapping not compitable for GA labels,(more than 4 taegets)')
        end if

        ! Reset fcReg to zero and apply at the right faces
        fcReg = 0
        indFc = (/ (i, i=1,grid%face%ntot) /)
        do i = 1, nl
            lb = lbls2(i)
            allocate(ind(count(fcLbl_loc == lb)))
            ind = pack(indFc,fcLbl_loc == lb )
            if (lb .lt. 1 .or. lb .gt. size(options%fcRegmappingGA)) then
                print *, 'Label out of range: ', lb
                call gdErrorHandler('PostProcessGA: label out of range for fcRegmappingGA, check facelabelmapping')
            end if
            fcReg(ind) = options%fcRegmappingGA(lb)
            deallocate(ind)
        end do

        ! Self-check if faces with fcReg label can be chained together
        do i = 1, size(options%fcRegmappingGA)

            if (options%fcRegmappingGA(i) /= 0) then
                nflbl = count(fcReg == options%fcRegmappingGA(i))
                allocate(fcs(nflbl))
                fcs = pack(indFc,fcReg == options%fcRegmappingGA(i))

                if (size(fcs) /= 0) then

                    ! Build chaines of faces
                    call grid%face%ChainFaces(fcs, f_ord, nf)
                    
                    ! One fcReg number should only have one chain
                    if (size(nf) > 1) then

                        call gdErrorHandler("PostprocessGA: more than one chain of faces detect with the fcReg value")

                    end if

                end if

                ! House keeping
                deallocate(fcs)


            end if
        end do

        ! Save
        call grid%face%reg%SetAllElementsArray(fcReg)

        ! Visualize end grid
        call grid%WriteData('grid_after_GA')


    end subroutine

    subroutine PostProcessingGridInformation(grid,magneticField,options)

        ! Description
        !============
        ! Post-process grid information to receive correct fields for
        ! traduit.out.b2us-file. This includes:
        ! - psi, bx, by, ffbz at the vertex locations
        ! - x, y coordinates of cells
        ! - bp, bt according to b2ag

        ! Notes
        !======
        ! Note 1: as implemented now, the grid.psi values may be inconsistent with
        ! the new grid coordinates, as first grid.psi is computed based on the old
        ! cell coordinates.  A similar issue is with the t0 term when recalculating 
        ! bp and bt. Needs revision        

        ! Declare variables
        !==================
        ! Arguments
        type(GridUDT), intent(inout)        :: grid
        type(MagneticFieldUDT), intent(in)  :: magneticField
        type(GAoptionsUDT), intent(inout)      :: options 

        ! Auxiliary
        logical :: trias_present, pents_present
        integer(I8):: vxs(1:100), vxsfd(grid%vert%ntot), s, nv, ic, i, ifs
        real(R8) :: bp(grid%cell%ntot), bt(grid%cell%ntot), &
            bpvx(grid%vert%ntot), fsPsi(grid%data%fluxdata%nFs), &
            r
        character(:), allocatable :: ctype

        ! Initialize
        associate(&
            c => grid%cell, &
            f => grid%face, &
            v => grid%vert &
            )

        ! Calculate additional connectivity
        call ComputeGridInterconnections(grid)

        ! Recalculate magnetic field
        call magneticField%interp%Evaluate(v%x,v%y,0,0,v%psi)
        call magneticField%interp%Evaluate(v%x,v%y,1,0,v%bx)
        call magneticField%interp%Evaluate(v%x,v%y,0,1,v%by)

        ! Recalculate ffbz (constant)
        v%ffbz = v%ffbz(1)
        if (magneticField%RBtor /= 0) then
            v%ffbz = 2*pi_R8*magneticField%RBtor
        end if 

        ! Cells
        !======
        ! Recalculate cell centers and bp and bt
        bp = 0
        bt = 0
        do ic = 1, c%ntot
            vxs = GetCellVert(c, ic)
            nv = c%vertP(ic,2)

            ! Compute cell centers as average of vertex coordinates
            c%x(ic) = sum(v%x(vxs(1:nv)))/real(nv, kind=R8)
            c%y(ic) = sum(v%y(vxs(1:nv)))/real(nv, kind=R8)

            ! Recompute bp and bt according to b2agbb (in b2ag (SOLPS-ITER))
            do i = 1, nv
                bpvx(i) = sqrt( v%bx(i)**2 + v%by(i)**2 )
            end do
            r = nv*2*pi_R8*c%x(ic)    ! For axissymmetry around Z = 0, x = 0 ! TODO
            bp(ic) = -sum(bpvx(1:nv))/r
            bt(ic) = sum(v%ffbz(vxs(1:nv)))/r
            bpvx = 0
        end do

        ! Save
        c%bp = bp
        c%bt = bt


        ! Recalculate psi at cell centers
        call magneticField%interp%Evaluate(c%x,c%y,0,0,c%psi)

        ! Get number of guard cells
        c%ngc = count(f%label /= 0)

        ! Ordening should be fine - done in PostProcessGA

        ! Build flux tube data
        !=====================
        trias_present = any(c%vertP(:,2) == 3)
        pents_present = any(c%vertP(:,2) == 5)

        if (.not.options%rem_small_trias) then
        
            print *, "Warning: Postprocessing: BuildFluxTube can not handle the" // & 
            & "presence of mini triangles as the face%aligned array is not correct"

        else 

            if (.not.trias_present .and. .not.pents_present) then

                ! For grids that can be stored as structured grid
                ctype = 'all'
                call BuildFluxTubeData(grid,options,magneticField,ctype)

            elseif (trias_present .and. .not. pents_present) then

                ! For grids with triangles (so fully unstructured) but without pentagons
                ctype = 'all_us'
                call BuildFluxTubeData(grid,options,magneticField,ctype); 

            elseif (pents_present) then

                ! For unstructured grids it is sufficient to give the closed flux tubes
                ! in the core and the first tube in the SOL of cells neighbouring the core. 
                ctype = 'closed'
                call BuildFluxTubeData(grid,options,magneticField,ctype);             

            end if
        end if 

        ! Recalculate fsPsi
        fsPsi = 0
        vxs = 0
        do ifs = 1, grid%data%fluxdata%nFs
            s = grid%data%fluxdata%fluxsurfacevertsP(ifs,1)
            nv = grid%data%fluxdata%fluxsurfacevertsP(ifs,2)
            vxsfd(1:nv) = grid%data%fluxdata%fluxsurfaceverts(s:s+nv-1)
            fsPsi(ifs) = sum(v%psi(vxsfd(1:nv))) / real(nv, kind=R8)
        end do

        ! Save
        grid%data%fluxdata%fluxsurfacepsi = fsPsi

        ! Determine OMP and IMP
        call DetermineMPs(grid, options)
        

        end associate

    end subroutine

    subroutine CarryOverOptions(goatoptions, gaoptions)

        ! Description
        !============
        ! Carry over some options from goatoptions to gaoptions

        ! Declare variables
        !==================
        ! Arguments
        type(GoatoptionsUDT), intent(in)    :: goatoptions
        type(GAoptionsUDT), intent(inout)   :: gaoptions
        
        gaoptions%vesselmode            = goatoptions%vesselmode 
        gaoptions%slab                  = goatoptions%slab
        gaoptions%debug                 = goatoptions%debug 
        gaoptions%facelabelmappingGG    = goatoptions%facelabelmappingGG
        gaoptions%facelabelmappingGA    = goatoptions%facelabelmappingGA
        gaoptions%facelabelmappingGD    = goatoptions%facelabelmappingGD 
        gaoptions%facelabelsubfrom      = goatoptions%facelabelsubfrom
        gaoptions%facelabelsubto        = goatoptions%facelabelsubto
        gaoptions%OMP_r                 = goatoptions%OMP_r
        gaoptions%OMP_z                 = goatoptions%OMP_z
        gaoptions%IMP_r                 = goatoptions%IMP_r
        gaoptions%IMP_z                 = goatoptions%IMP_z
        gaoptions%readstate             = goatoptions%readstate
        gaoptions%readstatemeth         = goatoptions%readstatemeth

    end subroutine

    subroutine CheckGAoptions(options)

        ! Description
        !============
        ! Check the consistency of the inputted options

        ! Declare variables
        !==================
        ! Arguments
        type(GAoptionsUDT), intent(inout) :: options

        ! Auxiliary
        integer(I8) :: nl, i

        ! BLG, first remove small triangles
        if (options%BLG) &
            options%rem_small_trias = .true.
        
        ! Pol flux
        if (options%rad_type == 3 &
            .and. options%splitting .and. .not.options%dist_function) then
            options%dist_function = .true.
            print *, 'Using pol_flux method for radial splitting while GAoptions.' // &
                & 'dist_function is off. Setting this to 1.'
            print *, 'options%dist_function: T'
        end if

        ! Make sure split and merge arrays are the same size
        nl = size(options%merging_array)
        if ( nl /= size(options%splitting_array) &
            .or. nl /= size(options%splittype_array) &
            .or. nl /= size(options%n_split_array) &
            .or. nl /= size(options%rad_type_array) &
            .or. nl /= size(options%pol_type_array) &
            .or. nl /= size(options%merge_crit_array) &
            .or. nl /= size(options%n_merge_array)) then
                call gdErrorHandler('CheckGAoptions: ga.splitting, ga.merging, ' // &
                    'ga.splittype, ga.n_split, ga.rad_type, ga.pol_type, ' // &
                    'ga.merge_crit, and ga.n_merge must have equal sizes')
        end if

        ! Enabled splitting passes must name a split direction explicitly.
        ! Other values would silently fall back to the default radial
        ! direction and bypass the checks that depend on the split type.
        do i = 1, nl
            if (options%splitting_array(i) == 1 .and. &
                options%splittype_array(i) /= 1 .and. &
                options%splittype_array(i) /= 2) then
                print *, 'CheckGAoptions: pass ', i, ' has ga.splittype = ', &
                    options%splittype_array(i)
                call gdErrorHandler('CheckGAoptions: ga.splittype must be ' // &
                    '1 (radial) or 2 (poloidal) for every enabled splitting pass')
            end if
        end do

        ! Aposteriori
        if (options%meth == 'aposteriori' .and. .not.options%readstate) then
            call gdErrorHandler('CheckGAoptions: to do aposteriori GA, a state need to read in')
        end if


    end subroutine

end module