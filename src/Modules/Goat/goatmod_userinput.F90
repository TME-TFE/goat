!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains all options and structures to process user input.
! It uses the modules mod_readwrite and mod_inputfileparser to process
! the input provided in the input file. 

module goatmod_userinput

    ! Initialize
    !============
    ! Load modules
    use mod_precision
    use mod_readwrite
    use mod_inputfileparser
    use mod_global_environment, only: solps, solps_inputfilepath, &
        solps_writefilepath, solps_gridfilepath, solps_magneticfieldfilepath, &
        solps_statefilepath, solps_structurefilepath, solps_outputfilepath

    ! The usual
    implicit none
    save
    public 

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                               GOAT                               !
    !------------------------------------------------------------------!

    ! General goat options
    type, extends(OptionsUDT) :: GoatoptionsUDT

        ! Structure containing the options for goat. The following 
        ! fields are present:

        ! General fields:
        ! - debug: general debug plots on the goat level
        ! - meth: goat running method. Currently, only 'GD' is supported
        ! - gridreadtype: type of grid input file to read. Can be 
        ! 'traduitb2us' or b2fgmtry (also in unstructured variant)
        ! - magneticfieldreadtype: type of magnetic field input to read.
        ! Can be 'rzpsi' or 'equ' for rzpsi.dat or equ.dat files.
        ! - filepath: path towards the file where options are defined
        ! - gdfilepath: path towards the file where options for grid 
        ! deformation are defined
        ! - readstate: logical to read state
        ! - readstatemeth: method to read state, can be 'b2fstate' or 'b2fplasmf'

        ! Input filenames:
        ! - gridfilepath: file path to file with grid data (e.g. traduit.out.b2us file)
        ! - structurefilepath: structure.dat file to read
        ! - magneticfieldfilepath: magnetic field file to read
        ! - statefilepath: state file to read
        ! - writefilepath: path where to write output traduit file

        ! Output options
        ! - write_final: write final output
        ! - write_traduitb2us: write unstructured traduit file
        ! - write_b2agdat:  write final b2ag.dat file for use in b2ag
        ! - write_OMPdata: write OMP data in traduit file
        ! - write_topologicaldata:   write out X-, O-, and strike point
        ! information, and other topological data. 

        ! Case identification options
        ! - vesselmode: set to true if the case is a vessel mode grid
        ! - slab: true if slab grid
        ! - artificial_slab: true if artificial slab

        ! Face label mappings
        ! - facelabelmappingGG: labels as defined in grid
        ! generator 
        ! - facelabelmappingGD: corresponding labels for GD (so 
        ! first GG label is mapped to first GD label here)
        ! - facelabelmappingGA: corresponding labels for GA (so 
        ! first GG label is mapped to first GA label here)
        ! - facelabelsubfrom: substitution of this face label ... 
        ! - facelabelsubto: ... to this face label in GA to GD 
        ! interface

        ! Structure options
        ! - TP: structure numbers that are target plates
        ! - TPind: indices of target plates
        ! - exclude: indices of structures in structure.dat that should 
        ! be excluded when generating the bounding polygon

        ! OMP and IMP definition
        ! - OMP_r, OMP_z: R, Z coordinates that define the outer mid 
        ! plane line segment
        ! - IMP_r, IMP_z: R, Z coordinates that define the inner mid 
        ! plane line segment

        ! Topological data
        

        ! General
        logical                     :: debug     
        character(:), allocatable   :: meth 
        character(:), allocatable   :: magneticfieldreadtype
        character(:), allocatable   :: filepath
        character(:), allocatable   :: gdinputfilepath
        logical                     :: readstate
        character(:), allocatable   :: readstatemeth


        ! Specify input filenames
        character(:), allocatable   :: gridfilepath
        character(:), allocatable   :: structurefilepath
        character(:), allocatable   :: magneticfieldfilepath
        character(:), allocatable   :: writefilepath

        ! Output options
        logical                     :: write_final 
        logical                     :: write_traduitb2us
        logical                     :: write_b2agdat
        logical                     :: write_OMPdata
        logical                     :: write_topologicaldata

        ! Case identification options
        logical                     :: vesselmode 
        logical                     :: slab 
        logical                     :: artificial_slab

        ! Face label mappings
        integer(I8), allocatable    :: facelabelmappingGG(:)
        integer(I8), allocatable    :: facelabelmappingGD(:) 
        integer(I8), allocatable    :: facelabelmappingGA(:)        
        integer(I8), allocatable    :: facelabelsubfrom(:) 
        integer(I8), allocatable    :: facelabelsubto(:) 
 

        ! Structure options
        integer(I8), allocatable    :: TP(:)
        integer(I8), allocatable    :: TPind(:) 
        integer(I8), allocatable    :: exclude(:)

        ! OMP and IMP
        real(R8), allocatable       :: OMP_r(:), OMP_z(:), IMP_r(:), &
            IMP_z(:)        

    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGoatOptions
        procedure   :: SetDefaults      => SetDefaultGoatOptions
        

    end type

    !------------------------------------------------------------------!
    !                          GRID ADAPTATION                         !
    !------------------------------------------------------------------! 
    type, extends(OptionsUDT) :: GAoptionsUDT   

        ! Structure containing the options for the grid adaptation part 
        ! of GOAT. The following fields are present: 

        ! General adaptation options
        ! - plt: enable general plotting
        ! - plt_qm: enable plotting quality metrics
        ! - meth: used method, can be 'simple' to adapt the grid based in 
        !         on grid metric, or can be 'aposteriori' which is using
        !         plasma state information

        ! Facelabel mapping
        ! - facelabelmappingGG: list of face labels of the grid generator
        ! - facelabelmappingGA: mapping of grid generator label to labels
        !                       required for grid adaptation
        ! - facelabelmappingGD: mapping of grid generator label to labels
        !                       required for grid deformation
        ! - facelabelsubfrom: face label of grid generator which should be 
        !                       substituted
        ! - facelabelsubto: substitution labels

        ! Operation options    
        ! - rem_small_trias: remove small triangles
        ! - cut_off_pol: threshold for small triangles poloidal height 
        !                   with respect to the poloidal neighbor
        ! - cut_off_surf: threshold for small triangles on the surface
        !                   area with respect to the surface area of its
        !                   neighbors
        ! - stacked_trias: apply transformation to stacked triangles
        ! - stacked_trias_checkAR: switch on selecting triangles based on 
        !                           its aspect ratio
        ! - stacked_trias_maxAR: maximal aspect ratio of a triangle
        !                        in the stacked triangle transformation
        ! - merge_stacked_trias: switch to apply merging stacked triangles
        !                        which are too skewed
        ! - merge_stacked_trias_angle_threshold: threshold of stacked triangle
        !                       too merge it (in degrees)
        ! - merge_trap_into_stacked: switch to allow merging a trapezoid
        !                            into a group of stacked triangles
        ! - stacked_to_cutcell: apply transformation from stacked triangles
        !                       to cutcell
        ! - stacked_to_cutcell: switch to transform to cutcell based on face
        !                       length of the aligned faces of the triangles
        ! - split_shaved_off_tube: apply radial splitting a concavely shaved
        !                           off fluxtube at the outer boundary
        ! - splitting: apply splitting of cells
        ! - merging: apply merging of cells
        ! - pents_to_tria: transfrom all remaining pentagon into triangles
        ! - h_rad_threshold: threshold of radial width of cell during merging,
        !                    only smaller cells are considered
        ! - h_rad_core_threshold: idem as above but only for cells in the core
        ! - BLG: add a boundary layer at the main targets of a grid
        ! - BLG_n_layers: determine the number of layers of the added 
        !                   boundary layer
        ! - BLG_rescaling_factor: size of the boundary layer wrt the 
        !                           upstream cells  
        ! - BLG_smoothing_factors: first number is the threshold ratio of 
        !              subsequent vertex displacement at which a rescaling need 
        !              to be done, the second number of the maximal ratio allowed (clip)
        ! - rem_stickout_trias: remove triangle that are only connected 
        !                       with the grid via one face
        ! - rem_trias_tube: remove flux tube with only two triangles
        ! - rem_outershell: remove or merge outer flux tubes with a triangle
        !                   at both ends
        ! - rem_tube_outershell_threshold: threshold for select a flux tube
        !                   to be merged with neighboring tube (h_rad of 
        !                   neighboring tube / h_rad of boundary tube)
        ! - outershell_handling: can be 'remove' or 'merge'
        ! - rem_stickout_quad: remove quad which are only connected to the grid
        !                       with one face
        ! - split_noalignedquads: splitting quads without aligned faces
        
        ! Splitting options        
        ! - no_pents: do not allow pentagons in the final grid
        ! - QTtype: method on how to handle radial splitting at inclined boundaries
        !           can be 'regular', i.e. just splitting a triangle radially, or
        !           can be 'pol-rad', i.e. continue poloidal splitting after
        !           the radial split an inclined boundary face
        ! - split_out: option to not split pentagons 
        ! - splittype: can be:
        !               1: for radial splitting, i.e. reducing radial width of cells
        !               2: for poloidal splitting, i.e. reducing the poloidal length of cells
        ! - n_split: number of allowed splitting operations      
        ! - typeT: method to split a triangle, can be 'stacked', i.e. triangle 
        !          is split into two trianlge, or 'cutcell', i.e. triangle is 
        !           split in a quad and a triangle
        ! - rad_type: criterium to apply radial splitting (see SelectSplitCell), can be:
        !             1: 'h_rad_psi': split cells with large psi width
        !             2: 'h_rad': split cells with large radial width
        !             3: 'pol_flux': split cells with large poloidal flux estimation based on 
        !                           pol_flux_est distance function
        !             4: 'farSOL': split highly inclined triangles in the farSOL
        !             5: 'farSOL_refinements': split cells with large psi width in the farSOL
        !             6: 'farSOLrefinement_targets': split cells with large psi width and 
        !                                         low boundary face inclination in the farSOL
        !             7: 'no_aligned_faces': split cells without aligned faces
        !             8: 'shaved_off_tubes': split cells in tube which is concavely shaved off
        ! - pol_type: criterium to aplly poloidal splitting (see SelectSplitCell) can be:
        !                1: 'trias_cvS': split triangles with largest surface area
        !                2: 'h_pol': split cells with largest poloidal length
        !                3: 'trias_farSOL': split triangles with largest surface area in the farSOL
        !                4: 'farSOLrefinement_hpol': split cells with largest poloidal length in the farSOL
        !                5: 'farSOLrefinement_targets': split cells in the farSOL with large
        !                                   poloidal length and low inclination at the boundary face
        ! - dist_function_threshold_split: threshold for the value the distance function 
        !               for selecting a cell to split. The value should be lower than 
        !               threshold.
        ! - dist_function_threshold_split_wall: idem for the distance function wrt the wall
        
        ! Merging options
        ! - no_hex: to not allow hexagonal cells in the final grid
        ! - merge_crit: criterium to select a face from which the two cell need to be merged
        !               can be (see SelectMergingFace) :
        !               1: 'tria_to_quad': merging two triangles to quadrilaterals
        !               2: 'min_area': merging cells which are smaller the mean surface area
        !               3: 'h_pol': merging cells with too small poloidal length
        !               4: 'bias': merging based on high bias between cells in poloidal direction
        !               5: 'pol_flux': merging cells with low poloidal flux estimation via 
        !                           similarly named distance function
        !               6: 'h_rad': merging cells with psi width smaller than h_rad_threshold
        !               7: 'h_rad_core': merging cells in the core radial width smaller than
        !                            h_rad_core_threshold
        !               8: 'bias_rad_farSOL': merging cells with large radial bias in farSOL
        !               9: 'bias_rad': merging cells with large radial bias
        !               10: 'skew_tria': merging triangles with low inclination and high aspect ratio
        ! - merge_h_pol_factor: factor to losen merge criterium based on poloidal length
        ! - n_merge: number of merging operations allowed
        ! - merge_bias_limit: thershold on bias to merge the cells. Merge is done when
        !                       real bias is larger than the threshold.
        ! - dist_function_threshold_merge: threshold for the value the distance function 
        !               for selecting a face to merge over. The value should be lower than 
        !               threshold.

        ! Pentagon options
        ! - no_pents_area_merge: use an area constraint on where pentagons are allowed during merging
        ! - no_pents_area_split: use an area constraint on where pentagons are allowed during splitting
        ! - no_pents_area_type: method to defined area where no pentagons are allowed. 
        !                       Can be 'coordinates' using no_pents_area_maxR etc. to define a box, 
        !                       or can be 'dist_function' to use a distance function. The value of 
        !                       the distance function should be lower the threshold (split or merge).
        ! - no_pents_area_maxR: option to recreate a box where not pentagons are allowed, 
        !                       idem for ..minR, ..maxZ, ..minZ
        ! - base_func: the basic function for the construction of the distance function, 
        !              can only be 'exp(-dist/d)'
        ! - dist_function: switch to construct distance functions
        ! - d_rescale: rescaling factor for characteristic length of distance function
        ! - d_rescale: idem as above but for wall distance function
        ! - dist_type: type of distance function can be 
        !              'target_single_null': use targets as reference 
        !              'target_to_vessel': use outer boundary of the mesh as reference
        !              'pol_flux_est': use separatrix as reference
        ! - dist_type_wall: type of distance function used for wall distance function, see above
        ! - d_char_type: characteristic length used to construct distance function, can be:
        !               'min_Xpoint_dist': minimal distance between targets and Xpoint
        !               'max_Xpoint_dist': maximal distance between targets and Xpoint

        ! Aposteriori adaptations
        ! - vxVol_style: style to interpolate from cell centers to vertices
        ! - apost_interpolation_meth: determines the interpolation method for interpolation, can be 'barycentric' or 'finite_element'
        ! - apost_interpolationC: level of continuity of interpolant
        ! - apost_interpolationM: order of interpolant
        ! - apost_meth: method to chose a cell to refine, can be 'grad'
        !   from read state and residual information to a newly adapted grid, or 
        !   can be 'res' to base the adaptation on residual information
        ! - apost_use_XX: flag to indicate to use a certain state field to decide
        !   which cells to split
        ! - apost_lambda_threshold: threshold for length scale per cell width, if length scale 
        !   per cell width is lower than the threshold, the cell is considered for refinement
        ! The rest of the options are carried over from goatoptions or are not changeable.

        ! General adaptation options
        logical                     :: plt
        logical                     :: plt_qm
        character(:), allocatable   :: meth


        ! Facelabel mapping
        integer(I8), allocatable    :: facelabelmappingGG(:)
        integer(I8), allocatable    :: facelabelmappingGA(:)
        integer(I8), allocatable    :: facelabelmappingGD(:)
        integer(I8), allocatable    :: facelabelsubfrom(:)
        integer(I8), allocatable    :: facelabelsubto(:)


        ! Operation options
        logical                     :: rem_small_trias
        real(R8)                    :: cut_off_pol
        real(R8)                    :: cut_off_surf

        logical                     :: stacked_trias
        logical                     :: stacked_trias_checkAR
        real(R8)                    :: stacked_trias_maxAR
        logical                     :: merge_stacked_trias
        real(R8)                    :: merge_stacked_trias_angle_threshold
        logical                     :: merge_trap_into_stacked

        logical                     :: stacked_to_cutcell
        logical                     :: stacked_to_cutcell_nonuniform

        logical                     :: split_shaved_off_tube


        logical                     :: splitting
        integer(I8), allocatable    :: splitting_array(:)
        logical                     :: merging
        integer(I8), allocatable    :: merging_array(:)
        logical                     :: pents_to_tria
        real(R8)                    :: h_rad_threshold
        real(R8)                    :: h_rad_core_threshold

        logical                     :: BLG
        integer(I8)                 :: BLG_n_layers
        real(R8)                    :: BLG_rescaling_factor
        real(R8), allocatable       :: BLG_smoothing_factors(:)

        logical                     :: rem_stickout_trias
        logical                     :: rem_trias_tube
        logical                     :: rem_outershell
        real(R8)                    :: rem_tube_outershell_threshold
        character(:), allocatable   :: outershell_handling
        logical                     :: rem_stickout_quad
        logical                     :: split_noalignedquads 
        
        ! Splitting options
        logical                     :: no_pents
        character(:), allocatable   :: QTtype
        logical                     :: split_out
        character(:), allocatable   :: splittype
        integer(I8), allocatable    :: splittype_array(:)
        integer(I8)                 :: n_split
        integer(I8), allocatable    :: n_split_array(:)
        character(:), allocatable   :: typeT
        integer(I8)                 :: rad_type
        integer(I8), allocatable    :: rad_type_array(:)
        integer(I8)                 :: pol_type
        integer(I8), allocatable    :: pol_type_array(:)
        real(R8)                    :: dist_function_threshold_split
        real(R8)                    :: dist_function_threshold_split_wall
        
        ! Merging options
        logical                     :: no_hex
        integer(I8)                 :: merge_crit
        integer(I8), allocatable    :: merge_crit_array(:)
        real(R8)                    :: merge_h_pol_factor 
        integer(I8)                 :: n_merge
        integer(I8), allocatable    :: n_merge_array(:)
        real(R8)                    :: merge_bias_limit
        real(R8)                    :: dist_function_threshold_merge    
        
        ! Pentagon options
        logical                     :: no_pents_area_merge
        logical                     :: no_pents_area_split
        character(:), allocatable   :: no_pents_area_type

        real(R8)                    :: no_pents_area_maxR
        real(R8)                    :: no_pents_area_minR
        real(R8)                    :: no_pents_area_maxZ
        real(R8)                    :: no_pents_area_minZ 
        
        ! Distance function options
        character(:), allocatable   :: base_func
        logical                     :: dist_function
        real(R8)                    :: d_rescale
        real(R8)                    :: d_rescale_wall
        character(:), allocatable   :: dist_type
        character(:), allocatable   :: dist_type_wall
        character(:), allocatable   :: d_char_type
        logical                     :: plt_dist_func  

        ! Aposteriori
        integer(I8)                 :: vxVol_style
        character(:), allocatable   :: apost_interpolation_meth
        integer(I8)                 :: apost_interpolationC
        integer(I8)                 :: apost_interpolationM
        character(:), allocatable   :: apost_meth
        real(R8)                    :: apost_lambda_threshold
        logical                     :: apost_use_na
        logical                     :: apost_use_ua
        logical                     :: apost_use_te
        logical                     :: apost_use_ti
        logical                     :: apost_use_tn
        logical                     :: apost_use_po
        logical                     :: apost_use_kt
        logical                     :: apost_use_zt
        logical                     :: apost_use_resco
        logical                     :: apost_use_resmo
        logical                     :: apost_use_resmt
        logical                     :: apost_use_reshe
        logical                     :: apost_use_reshi
        logical                     :: apost_use_reshn
        logical                     :: apost_use_respo
        logical                     :: apost_use_reszt
        logical                     :: apost_use_reskt

        
        
        ! Caried over from goatoptions
        ! Case identification
        logical                     :: debug        
        logical                     :: vesselmode
        logical                     :: slab
        logical                     :: readstate
        character(:), allocatable   :: readstatemeth    

        ! fcRegmappingGA
        integer(I8)                 :: fcRegmappingGA(1:7)

        real(R8), allocatable       :: OMP_r(:)
        real(R8), allocatable       :: OMP_z(:)
        real(R8), allocatable       :: IMP_r(:)
        real(R8), allocatable       :: IMP_z(:)

        ! Splitting
        logical                     :: XpointSplitting




    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGAOptions 
        procedure   :: SetDefaults      => SetDefaultGAoptions


    end type
   
    !------------------------------------------------------------------!
    !                          GRID DEFORMATION                        !
    !------------------------------------------------------------------!

    ! General grid deformation options
    type, extends(OptionsUDT) :: GDoptionsUDT

        ! Structure containing the options for the grid deformation part 
        ! of GOAT. The following fields are present:

        ! - runtype: type of run for grid deformation. currently, only 
        ! 'optimize' is available, though may be extended in the future
        ! - gridtype: only available options is plasma edge grid 
        ! ('plasma')
        ! - meth: only supported option is 'KKT', though augmented 
        ! lagrangian is also available. Method to solve optimization 
        ! problem when using 'optimize' as runtype

        ! - designoptionsfile: file to read design options
        ! - gridoptionsfile: similar to above, but for grid
        ! - magneticfieldoptionsfile: for MF
        ! - numparamsoptionsfile: for numerical parameters
        ! - environmentoptionsfile: for other environment stuff

        ! 
        ! - driver:     specify the driver to be used in the goat 
        !               program (see Goat.F90 for the options)
        ! - facelabelsubfrom:   labels in the original grid file to 
        !                       replace by other files
        ! - facelabelsubto:     replacement labels

        ! General deformation options
        character(:), allocatable   :: runtype 
        character(:), allocatable   :: gridtype 
        character(:), allocatable   :: meth 
        character(:), allocatable   :: filepath 

        ! Files to read other options
        character(:), allocatable   :: designoptionsfile 
        character(:), allocatable   :: gridoptionsfile 
        character(:), allocatable   :: magneticfieldoptionsfile 
        character(:), allocatable   :: numparamsoptionsfile
        character(:), allocatable   :: environmentoptionsfile

    contains

        ! Routines to manipulate the options
        procedure   :: Read             => ReadGDOptions
        procedure   :: SetDefaults      => SetDefaultGDOptions
        

    end type

    ! Options for the grid
    type, extends(OptionsUDT) :: GridOptionsUDT

        ! Structure containing options for grid manipulation (mainly
        ! reading and writing). The following fields are present:
        ! - type:       the type of grid that is considered. Only
        !                   'plasma' is currently available. 
        ! - readmeth:       the type of inputfile. Can be 'b2fgmtry' 
        !                   or 'traduit'. Both are unstructured!
        ! - filepath:       file path indicating where the grid file
        !                   data is stored

        ! Input filepaths and types
        character(:), allocatable   :: type 
        character(:), allocatable   :: readmeth 
        character(:), allocatable   :: filepath 

        ! Mappings
        integer(I8), allocatable    :: facelabelsubfrom(:), &
            facelabelsubto(:), facelabelmappingGG(:), &
            facelabelmappingGD(:), facelabelmappingGA(:)
        

    contains 

        procedure :: Read               => ReadGridOptions
        procedure :: SetDefaults        => SetDefaultGridOptions
        
    end type  

    ! Options for magnetic field
    type, extends(OptionsUDT) :: MagneticFieldOptionsUDT

        ! Structure containing options for reading in the magnetic
        ! field. The following fields are present:
        ! - readmeth:   method to read in the magnetic field data. Only
        !               'readrzpsi' and 'default' supported for now. 
        !               These use the readrzpsi routine
        ! - filepath:   path where the magnetic field data is stored
        ! - interpC:    desired continuity of the interpolant 
        ! - interpM     order of the interpolant describing the values 
        !               at the grid nodes 
        ! - interpmeth  interpolant method ('uniformgrid' for data given 
        !               on a uniform grid, 'centered' for default 
        !               non-uniform grids. 
        ! - RBtor       product of torodial magnetic field and major 
        !               radius (assumed constants here)
        ! - reinterpolate   switch to reinterpolate the magnetic field 
        !                   for a different resolution defined by resx, 
        !                   rexy

        
        character(:), allocatable   :: readmeth
        character(:), allocatable   :: filepath
        character(:), allocatable   :: interpmeth 

        integer(I8)                 :: interpC, interpM, resx, resy
        real(R8)                    :: RBtor
        logical                     :: reinterpolate 

    contains
    
        procedure :: Read           => ReadMagneticFieldOptions
        procedure :: SetDefaults    => SetDefaultMagneticFieldOptions
        
    end type  

    ! Options for the vessel
    type, extends(OptionsUDT) :: VesselOptionsUDT

        ! Fields:
        ! - readmeth:   'read_structure' for structure.dat files. Only
        !               method that is currently supported.
        ! - filepath:   path to the file to be read
        ! - refine:     set to true to do refinement of vessel (insert 
        !               more nodes)
        ! - maxdist:    maximum distance between two nodes. If larger,
        !               nodes will be added in between when refine == 1
        ! - minreffac   minimal refinement factor for each edge. If set 
        !               to zero, then no effect. Value+1 gives the number 
        !               of edges that each edge will be split minimally.
        ! - TP:         integer giving the indices of which part of the
        !               structure should be considered as target plates.
        !               This should not account for any structures being
        !               excluded using 'exclude' (see below)
        ! - TPind:      target plate enumerators 
        ! - exclude:    list of structures in the vessel to be excluded
        !               from the vessel 
        ! - shapemeth:  method used to represent (and possibly 
        ! approximate) the shape of the vessel. 
        ! - resx:       resolution in x-direction for shape rep
        ! - resy:       resolution in y-direction for shape rep
        ! - offsetfracx:    fractional offset to be taken from minimal 
        ! and maximal x-value of original vessel polygon
        ! - offsetfracy:    same as offsetfracx, but for y-direction
        ! - interpC:    only used for interpolant based methods. Order
        ! of interpolant
        ! - interpM:    only used for interpolant based methods. Order 
        ! used to approximate derivatives to construct interpolant 
        ! function. 
        ! - xrange, yrange: only used for interpolant based methods. 
        ! Sets range over which interpolant can be evaluated (range is 
        ! determined by max and min value in array, may be overwritten))

        character(:), allocatable       :: readmeth
        character(:), allocatable       :: filepath

        real(R8)                        :: maxdist
        logical                         :: refine
        integer(I8), allocatable        :: TP(:), TPind(:), exclude(:)
        integer(I8), allocatable        :: solpslimitertarget(:)
        integer(I8), allocatable        :: solpslowdivinnertarget(:)
        integer(I8), allocatable        :: solpslowdivoutertarget(:)
        integer(I8), allocatable        :: solpsupdivinnertarget(:)
        integer(I8), allocatable        :: solpsupdivoutertarget(:)

        ! Vessel representation options
        character(:), allocatable       :: shapemeth 
        integer(I8)                     :: resx, resy, interpC, interpM, &
            minreffac 
        real(R8)                        :: offsetfracx, offsetfracy
        real(R8), allocatable           :: xrange(:), yrange(:)

    contains 

        procedure :: Read           => ReadVesselOptions
        procedure :: SetDefaults    => SetDefaultVesselOptions
        
    end type 

    ! Options for the environment
    type, extends(OptionsUDT) :: EnvironmentOptionsUDT

        ! Structure used to keep all other possible necessary 
        ! structures that do not immediately fall under the grid, vessel
        ! or magnetic field. Currently empty. 
        ! - readstate: logical to read state
        ! - readstatemeth: method to read state, can be 'b2fstate' or 'b2fplasmf'

        character(:), allocatable   :: type
        character(:), allocatable   :: filepath
        character(:), allocatable   :: vesselfilepath
        logical                     :: readstate
        character(:), allocatable   :: readstatemeth
        character(:), allocatable   :: statefilepath

    contains 

        procedure :: Read           => ReadEnvironmentOptions
        procedure :: SetDefaults    => SetDefaultEnvironmentOptions
        
    end type

    !------------------------------------------------------------------!
    !                          GRID GENERATION                         !
    !------------------------------------------------------------------!

    ! Options for topological mesh generation
    type, extends(OptionsUDT) :: TopomeshOptionsUDT 

        ! Description
        !============
        ! Options for the topological mesh construction. The 
        ! following fields are set:
        ! - fresx, fresy    : resolution for 2D tracing of contours and
        !                   extrema of the field
        ! - fdonewton       : option to refine extrema with newton 
        !                   solver (may not always converge!)
        ! - dotpvesselbased : do tangency point determination purely 
        !                   based on current vessel polygons 
        ! - v(...)          : same options but for vessel 
        ! - ffieldtol       : tolerance on field value of extrema (if 
        !                   difference is below tolerance, two extrema 
        !                   are considered to have exactly the same field
        !                   value)
        ! - doadaptations       : general switch to apply or not apply
        !                       adaptations to the basic topological 
        !                       mesh (all things like core/PF boundaries, 
        !                       optional region removal, tube merging, ...)
        ! - addcoreboundaries   : adds additional core boundaries
        ! - addPFboundaries     : adds additional private flux like 
        !                       boundaries (at tangency points of which 
        !                       contour doesn't enter vessel)
        ! - removecoreregions   : remove the innermost core regions
        ! - coreboundariesfrac  : fraction in field value between core 
        !                       and connecting X-point (1: tangency point, 0: 
        !                       other PF boundary, e.g. separatrix)
        ! - PFboundariesfrac    : fraction in field value between outer 
        !                       PF boundary and tangency point
        ! - removewidegridregions: remove all regions that are not next    
        !                           to a separatrix 
        ! - removevesselregions : remove regions that are adjacent to 
        !                       some user specified vessel region
        ! - rvrvesselIDs        : IDs, as specified in the structure.dat
        !                       or vessel.dat file, to consider
        ! - rvrretain:          : retain tubes instead of deleting them, 
        !                       and delete all others
        ! - rvrfullycovered     : if true, a face is only considered if 
        !                       it is fully covered by one (or multiple) 
        !                       defined vessel IDs
        ! - rvrdocascade        : delete not only marked tubes, but also
        !                       tubes that connect (so cascade the deletion)
        !                       in either increasing or decreasing psi 
        !                       direction (next option)
        ! - rvrcascadedir       : 'upwards' for tubes with higher psi 
        !                       value, 'downwards' for lower, 'none' for 
        !                       no cascade. 
        ! - npmin, npmax, dl    : minimal and maximal number of points
        !                       of contours (when doing coarsening) and
        !                       desired uniform edge length
        ! - readexistingTM:     read in an existing topomesh file, 
        !                       for which the full path is defined in  
        !                       TMfilepath
        ! - readexistingtracers     read in existing field and vessel
        !                       tracers instead of constructing
        !                       new ones - only when restarting the 
        !                       topomesh from a previous one. field and
        !                       vessel filepaths should be fully specified
        !                       in TMfieldtracerfilepath and TMvesseltracerfilepath
        ! - writedebugoutput    write (many) more intermediate topomesh 
        !                       files for visual debuggging (default false)
        ! - mergetubes          apply general merging (should outperform 
        !                       mergetangencypointtubes and mergeavptubes)
        ! - mergetubemeth       method to merge tubes, either 'extensive' or 'local'. 
        !                       the extensive method considers several 
        !                       tubes at the same time and often leads 
        !                       to 'nicer' tubes for grid generation, but
        !                       may not be able to merge some tubes when 
        !                       when cycles are present in the tube graph.
        !                       the local method considers only two tubes
        !                       at the same time and can therefore merge    
        !                       more tube pairs, but may lead to uglier tubes. 
        ! - dpsimintubes        minimal psi difference for topomesh tubes (below
        !                       which we try to merge - set to 0 to ignore)
        ! - lradmintubes        minimal radial length for topomesh tubes
        !                       (below which we try to merge - set to 0 to ignore)
        ! - lpolsepmin          minimal separatrix face length in poloidal direction.
        !                       If lower, separatrix faces can be removed
        !                       during separatrix simplification step
        ! - mtallowseparatrix   allow merging over the separatrix for 
        !                       the general merge case
        ! - mtallowcore, mtallowpf  same as separatrix, but for core and 
        !                       pf contours
        ! - alignvesselparts    define certain vessel parts as aligned 
        !                       faces with a certain flux surface value 
        !                       and flux surface ID. Only certain boundary
        !                       faces will be considered for alignment 
        !                       (typically those near type 2 tangency 
        !                       points). This will only be done after
        !                       initial topological mesh construction 
        !                       as this triggers profound adaptation of     
        !                       of the topomesh.
        ! - avpminangle         minimum angle w.r.t. the magnetic field 
        !                       of boundary edges. If below, it will be
        !                       considered as potential aligned part 
        !                       (avp: aligned vessel parts). This angle
        !                       should be given in degrees!
        ! - (deprecated) avprefinevessel     switch to refine vessel boundaries, similar
        !                       to full vessel refinement (see vessel options)
        ! - (deprecated) avpmaxvesseldist    maximal vessel edge length
        ! - (deprecated) avpminreffac        minimal refinement factor for vessel edge refinement                    
        ! 

        integer(I8)             :: fresx, fresy, vresx, vresy, npmin, &
            npmax, avpminreffac
        integer(I8), allocatable, dimension(:)  :: rvrvesselIDs
        logical                 :: addcoreboundaries, removecoreregions, &
            fdonewton, vdonewton, removewidegridregions, addPFboundaries, &
            readexistingTM, removenoncoreregions, &
            doadaptations, dotpvesselbased, removevesselregions, rvrretain, &
            rvrdocascade, rvrfullycovered, alignvesselparts, avprefinevessel, &
            readexistingtracers, mergetubes, mtallowcore, simplifyseparatrix, &
            mtallowpf, mtallowseparatrix, writedebugoutput, mtoldstyle
        real(R8)                :: coreboundariesfrac, ffieldtol, dl, &
            PFboundariesfrac, avpminangle, avpmaxvesseldist, &
                dpsimintubes, lradmintubes, lpolsepmin
        character(:), allocatable   :: TMfilepath, rvrcascadedir, &
            TMfieldtracerfilepath, TMvesseltracerfilepath, mergetubemeth
    contains 

        procedure :: Read           => ReadTopomeshOptions
        procedure :: SetDefaults    => SetDefaultTopomeshOptions

    end type

    ! Options for grid generation
    type, extends(OptionsUDT) :: GGOptionsUDT 

        ! Description
        !============
        ! Options for the grid generator, including vertex distribution 
        ! options (poloidal and radial), flux surface removal options, 
        ! boundary triangle removal options, and small face removal 
        ! options. The following fields are present: 

        ! General grid generation
        ! - verbosity       the higher, the more information is 
        !                   printed out (default: 1, 0 suppresses all)
        ! - ggmethod        method how to construct the grid. Can be 
        !                   'independent' (treating each flux surface
        !                    independently) or 'orthogonal' (
        !                   yields largely orthogonal grid, but needs 
        !                   dependency between mesh cells)
        ! - cellconstructionmethod: method how to determine grid cells
        !                   starting from given vertex distribution. 
        !                   'quads_triangles' is recommended one
        ! - TMcellgriddingorder:    order in which to grid the 
        !                   topological mesh cells. Can be 'independent'
        !                   (then it is  pretty much random) or 
        !                   'sequential' -> that one is the recommended
        !                   one. Here, initial distributions are 
        !                   propagated through
        ! 
        ! Vertex distribution, poloidal direction:
        ! - vdptype:        'uniform' for uniform distribution, 
        !                   'densitybased' for non-uniform distribution 
        !                   based on predefined node density distribution
        !                   function 
        ! - vdpdfacelength: length of face for 'uniform' distribution 
        ! - vdpddecaylength:    decay length from vessel for 
        !                       'densitybased' option
        ! - vdpdx, y, d, val:   parameters to add attraction points for
        !                       densitybased option: x, y are 
        !                       coordinates, d is a decay length, and 
        !                       val is the desired density at the point
        ! - vdpddensityatvessel:    desired density at the vessel boundary
        ! - vdpddensityatinf:       desired density far from the vessel 
        !                           boundary 

        ! Vertex distribution, radial direction:
        ! largely the same as poloidal one (but then vdr instead of vdp)
        ! but some differences
        ! - vdrtype:        (same as vdptype)
        ! - vdrdoxp:            include x-points? (for refinement in core region)
        ! - vdrdosp:        include strike points? (for refinement near targets)
        ! - vdrdfieldwidth: desired field with for uniform distribution
        ! - vdrddecaylength:    decay length parameter 
        ! - vdrddensityatxp:    desired density at separatrix
        ! - vdrddensityatinf:   density far from separatrix
 
        ! Options for extending flux tubes with vessel parts (so-called 
        ! 'cut cells')
        ! - extendtptubes:          extend tubes at the tangency point
        !                           side to include (part of) the 
        !                           vessel from both side of type 1
        !                           tangency points
        ! - extendvesseltubes:      extend any tube at the vessel edges
        !                           that complies to the marking criterion
        ! - evtmaxvessellength:     maximum L2-based vessel segment 
        !                           length before tube is extended
        ! - evtnoBL                 if true, don't apply any boundary 
        !                           layer at extended tubes (typically 
        !                           not desired. Only has effect if 
        !                           boundary layers are applied of course)

        ! Options for flux surface removal
        ! - removefluxsurfaces:     switch to remove or not
        ! - remfspsitol:            absolute tolerance for difference in 
        !                           flux value between flux surfaces 
        !                           (if below, one is removed)
        ! - remfspsirattol:         maximal ratio between flux surface 
        !                           difference (and inverse), if above/
        !                           below: removed

        ! Options for boundary triangle removal
        ! - removenarrowboundarytriangles:  switch
        ! - rembndtùriacriterion:    criterion for removal, typically 
        !                           'angle' (if too small, remove triangle)
        ! - rmbndtriaminangle:      minimal angle [rad] (too small - removed)
        !                           input is in degrees!
        
        ! Options for small face removal
        ! - removefaces             switch 
        ! - remfacescriterion       criterion to remove faces. Typically
        !                           'facelength_radial_bnd' to remove
        !                           too narrow radial faces at the 
        !                           boundary 
        ! - remfacesminlength       minimal length of these faces  [m] 

        ! Options for refinement
        ! - refmeth:                method for refinement. Can be 'no' 
        !                           (no additional refinement), 
        !                           'lengthbased' (ref based on min and
        !                           max length distributions) 
        ! - reflengthtype:          which type of length to consider. 
        !                           'euler' is classical eulerian length 
        !                           (L2 norm), 'radial' is projected in 
        !                           radial direction (and absolute value taken)
        
        ! (poloidal) Refinement options for lengthbased option (this is currently 
        ! based on exponential decay functions defined in points):
        ! - refLBlmininf    minimal length at infinity
        ! - refLBLmaxinf    maximal length at infinity
        ! - refLBdoxp:      refine near x-points (x-points are added)
        ! - refLBLminxp     minimal length on x-point
        ! - refLBLmaxxp     maximal length on x-point
        ! - refLBdecaylengthxp  decaylength on x-point (larger - wider influence)
        ! - refLBdovessel   refine near (specified) vessel vertices
        ! - refLBstructureIDs   structure IDs of vessel structures to include
        ! - refLBvertIDs:   vessel vertex IDs of vessel vertices to include, 
        !                   as specified in structure.dat file (some 
        !                   may not be included if deleted during vessel 
        !                   polygon construction)
        ! - refLBLminstructure  minimal length on structure i (etc, 
        !                   similar for vertices)

        ! - refBLdotarget   do BL refinement at targets
        ! - refBLdovessel   do BL refinement at far vessel boundaries
        ! - refBLdostructure    do BL refinement on user-specified structures
        ! - refBLnctarget   number of desired boundary layer cells at 
        !                   the target (similar for vessel and structure)
        ! - refBLdltarget   desired lengths for these cells 
        ! - refdlBLlengthbased  desired length is specified in classic 
        !                       euler length or not (if not, lengthtype 
        !                       is taken)
        ! - refBLstructureID    specific for structure based BL: the 
        !                       structure labels (as specified in the   
        !                       input structure file) at which refinement
        !                       should be done 

        
        ! - refBLncstructure    number

        ! (Radial) refinement options for lengthbased refiner:
        !   mostly the same refinement options as the poloidal direction,
        !   but names now have 'rad' in front. 
        ! - radrefLBlmininf    minimal length at infinity (in [m]!)
        ! - radrefLBLmaxinf    maximal length at infinity
        ! - radrefLBdosp:      refine near strike points (x-points are added)
        ! - radrefLBLminsp     minimal length on strike point
        ! - radrefLBLmaxsp     maximal length on strike point
        ! - radrefLBdecaylengthsp  decaylength on strike point (larger - wider influence)
        
        ! - radrefBLdosp    do BL refinement at strike points
        ! - radrefBLncsp   number of desired boundary layer cells at 
        !                   the strike point 
        ! - radrefBLdlsp   desired lengths for these cells

        ! Cell distribution options
        ! - legalcellstyle      option to control how to check if cells
        !                       are legal. if 'no', then no checks are
        !                       made. 'old' is an old, deprecated way 
        !                       that does not always work but catches 
        !                       most of the issues. If neither of the
        !                       previous options are taken, it defaults
        !                       to a graph-based method that is likely
        !                       more expensive, but does capture all 
        !                       possible cell overlap cases (normally
        !                       speaking)

        ! Label translation options:
        !   - structurebasedlabels:     base labels on structure IDs 
        !   - forceSOLPStopology:       force certain topology for region translation etc
        !   - SOLPStopology:            desired topology ('general', 'linear',
        !                               'limiter', 'single_null', 'double_null')

        ! Diagnostics
        ! - dogriddiagnostics   run grid diagnostics. Will be time consuming!

        
        logical                     :: removefluxsurfaces, &
            removenarrowboundarytriangles, removefaces, refLBdoxp, &
            refLBdovessel, vdpdincludexp, coarsencontours, refBLdotarget, &
            refBLdovessel, readexistingrefdata, radrefBLdosp, radrefLBdosp, &
            extendtptubes, extendvesseltubes, refdlBLlengthbased, &
            radrefdlBLlengthbased, vdrdoxp, structurebasedlabels, &
            dogriddiagnostics, evtnoBL, refBLdostructure, &
            forceSOLPStopology, vdrdosp, &
            solpstargetrefinements
        integer(I8)                 :: gcresx, gcresy, &
            verbosity, orthtracernsteps, refBLnctarget, refBLncvessel, &
            radrefBLncsp, refBLncstructure
        integer(I8), allocatable, dimension(:)  :: refLBstructureIDs, &
            refLBvertIDs, refBLstructureID
        real(R8)                    :: vdpdfacelength, vdpddecaylengthplf, &
            vdpddecaylengthxp, vdpddensityatvessel, vdpddensityatxp, &
            vdpddensityatinf, vdrdfieldwidth, vdrdfacelength, &
            vdrddecaylengthxp, vdrddensityatxp, vdrddensityatinf, &
            vdrddecaylengthsp, vdrddensityatsp, &
            remfspsitol, remfspsirattol, rembndtriaminangle, &
            remfacesminlength, refLBLmininf, refLBLmaxinf, refLBLminxp, &
            refLBLmaxxp, refLBdecaylengthxp, orthtracerstep, &
            radrefLBLmininf, radrefLBLmaxinf, radrefLBLminsp, &
            radrefLBLmaxsp, radrefLBdecaylengthsp, evtmaxvessellength
        real(R8), allocatable, dimension(:)     :: vdpdx, vdpdy, vdpdd, &
            vdpdval, refLBLminstructure, refLBLminvert, refLBLmaxstructure, &
            refLBLmaxvert, refLBdecaylengthstructure, refLBdecaylengthvert, &
            refBLdltarget, refBLdlvessel, radrefBLdlsp, vdrdx, vdrdy, &
            vdrdd, vdrdval, refBLdlstructure
        character(:), allocatable   :: vdptype, vdpdtype, vdrtype, &
            vdrdtype, rembndtriacriterion, remfacescriterion, ggmethod, &
            cellconstructionmethod, TMcellgriddingorder, refmeth, vdpplftype, &
            refdatafile, radrefmeth, reflengthtype, radreflengthtype, &
            legalcellstyle, SOLPStopology
    contains 

        procedure :: Read           => ReadGGOptions
        procedure :: SetDefaults    => SetDefaultGGOptions

    end type


    contains 

    !==================================================================!
    !                                                                  !
    !                               ROUTINES                           !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                            Option setters                        !
    !------------------------------------------------------------------!

    ! Goat options routines
    subroutine SetDefaultGoatOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)       :: options        
        
        ! Input file
        if (solps) then 
            ! SOLPS defaults
            options%inputfilepath    = solps_inputfilepath
        else
            options%inputfilepath    = './GOAToptions.dat'
        end if 

        ! General
        options%debug           = .false. 
        options%meth            = 'GD'
        options%gdinputfilepath = './GOAToptions.dat'

        ! Specify input filenames
        if (solps) then 
            ! SOLPS defaults
            options%gridfilepath            = solps_gridfilepath
            options%structurefilepath       = solps_structurefilepath
            options%magneticfieldfilepath   = solps_magneticfieldfilepath
        else
            options%gridfilepath            = './traduit.out.b2us'
            options%structurefilepath       = './structure.dat'
            options%magneticfieldfilepath   = './rzpsi.dat'
        end if 
        

        ! Output options
        if (solps) then 
            ! SOLPS defaults
            options%writefilepath   = solps_writefilepath
        else
            options%writefilepath       = 'traduit.out.b2us_smoothed'
        end if 
        options%write_final         = .true. 
        options%write_traduitb2us   = .true.
        options%write_b2agdat       = .true. 
        options%write_OMPdata       = .false. 
        options%write_topologicaldata = .false.

        ! Case identification options
        options%vesselmode          = .false. 
        options%slab                = .false.
        options%artificial_slab     = .false.
        
        ! Face label mappings
        allocate(options%facelabelmappingGG(0), &
            options%facelabelmappingGD(0), &
            options%facelabelmappingGA(0), &
            options%facelabelsubfrom(0), &
            options%facelabelsubto(0))
        
        ! Structure options
        allocate(options%TP(0), options%TPind(0), options%exclude(0))

        ! OMP and IMP
        allocate(options%OMP_r(1:2), options%OMP_z(1:2), &
            options%IMP_r(1:2), options%IMP_z(1:2))
        options%OMP_r(:)   = 0
        options%OMP_z(:)   = 0
        options%IMP_r(:)   = 0
        options%IMP_z(:)   = 0

    end subroutine

    ! Grid adaptations options routines
    subroutine SetDefaultGAoptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GAoptionsUDT)       :: options   
        
        ! General options
        options%plt         = .false.
        options%plt_qm      = .false.
        options%meth        = 'simple'

        ! Set fcReg mapping
        options%fcRegmappingGA = [0, 0, 0, 1, 4, 5, 8]

        ! Operation options
        options%rem_small_trias                     = .false.
        options%cut_off_pol                         = 0.3
        options%cut_off_surf                        = 0.05
    
        options%stacked_trias                       = .false.
        options%stacked_trias_checkAR               = .false.
        options%stacked_trias_maxAR                 = 10
        options%merge_stacked_trias                 = .false.
        options%merge_stacked_trias_angle_threshold = 5
        options%merge_trap_into_stacked             = .true.

        options%stacked_to_cutcell                  = .false.
        options%stacked_to_cutcell_nonuniform          = .true.

        options%split_shaved_off_tube               = .false.

        options%splitting                           = .false.
        options%splitting_array                     = [0]
        options%pents_to_tria                       = .false.
        options%merging                             = .false.
        options%merging_array                       = [0]
        options%h_rad_threshold                     = 0.01
        options%h_rad_core_threshold                = 0.04

        options%BLG                                 = .false.
        options%BLG_n_layers                        = 0
        options%BLG_rescaling_factor                = 2
        options%BLG_smoothing_factors               = [1.2, 1.5]

        options%rem_stickout_trias                  = .false.
        options%rem_trias_tube                      = .false.
        options%rem_outershell                      = .false.
        options%rem_tube_outershell_threshold       = 2
        options%outershell_handling                 = 'merge' 
        options%rem_stickout_quad                   = .false.
        options%split_noalignedquads                = .true. 
        
        ! Splitting options
        options%no_pents                            = .true.
        options%QTtype                              = 'regular'
        options%split_out                           = .false.
        options%splittype                           = 'rad'
        options%n_split                             = 20
        options%n_split_array                       = [20]
        options%typeT                               = 'cutcell'
        options%rad_type                            = 1
        options%pol_type                            = 1
        options%dist_function_threshold_split       = 0.9
        options%dist_function_threshold_split_wall  = 0.6  
        
        ! Merging options
        options%no_hex                              = .true.
        options%merge_crit                          = 4
        options%merge_crit_array                    = [4]
        options%merge_h_pol_factor                  = 1
        options%n_merge                             = 20
        options%n_merge_array                       = [20]
        options%merge_bias_limit                    = 5
        options%dist_function_threshold_merge       = 0.6   
        
        ! Pentagon options
        options%no_pents_area_merge                 = .false.
        options%no_pents_area_split                 = .false.
        options%no_pents_area_type                  = 'dist_function'
        
        options%no_pents_area_maxR                  = 2.5
        options%no_pents_area_minR                  = 1
        options%no_pents_area_maxZ                  = -0.9
        options%no_pents_area_minZ                  = -2
        
        ! Distance function options
        options%base_func                           = 'exp(-dist/d)'
        options%dist_function                       = .true.
        options%d_rescale                           = 0.5
        options%d_rescale_wall                      = 0.5
        options%dist_type                           = 'target_single_null'
        options%dist_type_wall                      = 'target_to_vessel'
        options%d_char_type                         = 'min_Xpoint_dist'
        options%plt_dist_func                       = .false.    

        ! Aposteriori
        options%vxVol_style                         = 2
        options%apost_interpolation_meth            = 'barycentric'
        options%apost_interpolationC                = 0
        options%apost_interpolationM                = 1
        options%apost_meth                          = 'grad'
        options%apost_lambda_threshold              = 3
        options%apost_use_na                        = .false.
        options%apost_use_ua                        = .false.
        options%apost_use_te                        = .false.
        options%apost_use_ti                        = .false.
        options%apost_use_tn                        = .false.
        options%apost_use_po                        = .false.
        options%apost_use_kt                        = .false.
        options%apost_use_zt                        = .false.
        options%apost_use_resco                     = .false.
        options%apost_use_resmo                     = .false.
        options%apost_use_resmt                     = .false.
        options%apost_use_reshe                     = .false.
        options%apost_use_reshi                     = .false.
        options%apost_use_reshn                     = .false.
        options%apost_use_respo                     = .false.
        options%apost_use_reszt                     = .false.
        options%apost_use_reskt                     = .false.

    end subroutine


    ! Grid deformation options routines
    subroutine SetDefaultGDOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(GDoptionsUDT)       :: options        
        
        ! Input file
        options%filepath    = './GOAToptions.dat'

        ! General
        options%runtype     = 'optimize'
        options%gridtype    = 'plasma'
        options%meth        = 'KKT'

        ! Files to read other options
        options%designoptionsfile           = './GOAToptions.dat'
        options%gridoptionsfile             = './GOAToptions.dat'
        options%magneticfieldoptionsfile    = './GOAToptions.dat'
        options%numparamsoptionsfile        = './GOAToptions.dat'
        options%environmentoptionsfile      = './GOAToptions.dat'
        
    end subroutine

    ! Grid option routines
    subroutine SetDefaultGridOptions(options)

        ! Description
        !============
        ! Set the options for the grid

        ! Declaration
        class(GridOptionsUDT)        ::  options

        ! Default options
        options%type                = 'plasma'
        options%readmeth            = 'b2fgmtry'

        ! Default grid location
        options%filepath            = 'traduit.out.b2us'

        ! Default mappings
        if (.not.allocated(options%facelabelsubfrom)) &
            allocate(options%facelabelsubfrom(0), options%facelabelsubto(0))
        if (.not.allocated(options%facelabelmappingGG)) &
            allocate(options%facelabelmappingGG(1:8), options%facelabelmappingGD(1:8), &
        options%facelabelmappingGA(1:8))
        options%facelabelmappingGG = [-13, -34, -23, -24, -21, -42, -43, -44]
        options%facelabelmappingGD = [1, 2, 3,   3,   4,   5,   5,   5]
        options%facelabelmappingGA = [4, 5, 3,   3,   2,   3,   3,   3]        

    
    end subroutine

    ! Magnetic field
    subroutine SetDefaultMagneticFieldOptions(options)
        ! Description
        !============
        ! Set the options for the magnetic field. The following fields
        ! have to be set:

        ! Declaration
        class(MagneticFieldOptionsUDT)       :: options

        ! Default options   
        options%readmeth                = 'default'
        options%filepath                = './Examples/TCV/rzpsi_tcv.dat'

        options%interpmeth              = 'uniformgrid' 
        options%interpC                 = 3
        options%interpM                 = 6
        options%reinterpolate           = .false. 
        options%resx                    = 100 
        options%resy                    = 100

        options%RBtor                   = 0

    end subroutine

    ! Vessel
    subroutine SetDefaultVesselOptions(options)

        ! Description
        !============
        ! Set the options for the vessel.

        ! Declare variables
        !==================
        class(VesselOptionsUDT)          :: options

        ! Set options
        !============
        ! Input
        options%readmeth    = 'read_structure'
        options%filepath    = './structure.dat'

        ! Refinement options
        options%refine      = .false.
        options%maxdist     = 0.01
        options%minreffac   = 0_I8
 
        ! Target plates
        allocate(options%TP(2))
        allocate(options%TPind(2))
        allocate(options%exclude(0))
        allocate(options%xrange(0), options%yrange(0))
        options%TP      = [1, 2]
        options%TPind   = [1, 2]

        ! SOLPS limiter target segments (empty by default; required
        ! non-empty for SOLPS limiter grids)
        allocate(options%solpslimitertarget(0))

        ! SOLPS divertor target segments (empty by default; consumed
        ! by the SOLPS region assignment of divertor topologies)
        allocate(options%solpslowdivinnertarget(0))
        allocate(options%solpslowdivoutertarget(0))
        allocate(options%solpsupdivinnertarget(0))
        allocate(options%solpsupdivoutertarget(0))

        ! Vessel shape representation
        allocate(character(len('closedpolygon_smoothapproximation')) :: &
            options%shapemeth)
        options%shapemeth = 'closedpolygon_smoothapproximation'
        options%resx            = 100
        options%resy            = 100 
        options%offsetfracx     = 0.1
        options%offsetfracy     = 0.1
        options%interpC         = 3
        options%interpM         = 6
        
        

    end subroutine

    ! Environment
    subroutine SetDefaultEnvironmentOptions(options)

        ! Description
        !============
        ! Set the options to read in the environment.

        ! Declare variables
        !==================
        class(EnvironmentOptionsUDT)         :: options

        ! Set options
        !============
        options%type = 'vessel'
        options%vesselfilepath = options%inputfilepath

        ! Carry over some goatoptions to environmentoptions
        options%readstate       = .false.
        options%readstatemeth   = 'b2fplasmf'
        
        ! Specify input filename
        if (solps) then
            ! SOLPS default
            options%statefilepath           = solps_statefilepath
        else
            options%statefilepath           = './b2fplasmf'
        end if

        ! Set data

    end subroutine

    ! Topomesh
    subroutine SetDefaultTopomeshOptions(options)

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshOptionsUDT)   :: options

        ! Set defaults
        !=============
        ! I/O
        options%readexistingTM = .false. 
        options%readexistingtracers = .false. 
        options%TMfilepath = 'topomesh.dat'
        options%TMfieldtracerfilepath = './output/TMfieldtracer.dat'
        options%TMvesseltracerfilepath = './output/TMvesseltracer.dat'
        options%writedebugoutput = .false.
        
        ! Contouring (field)
        options%fresx = 100
        options%fresy = 100

        ! Contouring (vessel)
        options%vresx = 100
        options%vresy = 100

        ! Contouring (general)
        options%npmin = 10
        options%npmax = 10000
        options%dl = 1e-3

        ! Refining options for extrema (field)
        options%fdonewton = .true. 
        options%ffieldtol = 1e-8 
        
        ! Refining options for extrema (vessel)
        options%vdonewton = .true.

        ! Additional options
        options%dotpvesselbased             = .false.

        ! Adaptations
        options%doadaptations               = .true.
        options%addcoreboundaries           = .true. 
        options%coreboundariesfrac          = 0.2
        options%addPFboundaries             = .true. 
        options%PFboundariesfrac            = 0.2
        options%removecoreregions           = .true. 
        options%removewidegridregions       = .true. 
        options%removenoncoreregions        = .false.
        options%mergetubes                  = .false. 
        options%mergetubemeth               = 'extensive'
        options%mtallowcore                 = .false.
        options%mtallowpf                   = .false.
        options%mtallowseparatrix           = .false.
        options%mtoldstyle                  = .false. 
        options%simplifyseparatrix          = .false. 
        options%dpsimintubes                = 0.0_R8 ! zero to ignore
        options%lradmintubes                = 0.0_R8 ! zero to ignore
        options%lpolsepmin                  = 0.0_R8

        options%removevesselregions         = .false. 
        if (allocated(options%rvrvesselIDs)) deallocate(options%rvrvesselIDs)
        allocate(options%rvrvesselIDs(0))
        options%rvrfullycovered             = .false. 
        options%rvrretain                   = .false. 
        options%rvrdocascade                = .false. 
        options%rvrcascadedir               = 'none'

        options%alignvesselparts            = .false.
        options%avprefinevessel             = .false.  
        options%avpminangle                 = 0.0_R8
        options%avpmaxvesseldist            = 1e-2 ! [m]
        options%avpminreffac                = 0

    end subroutine 

    ! Grid generation
    subroutine SetDefaultGGOptions(options)

        ! Modules
        !========
        use mod_constants, only : pi_R8

        ! Declare variables
        !==================
        ! Arguments
        class(GGOptionsUDT)   :: options

        ! Set defaults
        !=============
        ! Contouring options in grid generator
        options%gcresx = 100 
        options%gcresy = 100
        options%coarsencontours = .false.

        ! Options for poloidal vertex distribution
        options%vdptype             = 'densitybased'
        options%vdpplftype          = 'target'
        options%vdpdfacelength      = 4e-2 
        options%vdpddecaylengthplf  = 0.05
        options%vdpddecaylengthxp   = 0.05
        options%vdpdincludexp       = .false.
        allocate(options%vdpdx(0), options%vdpdy(0), options%vdpdd(0), &
            options%vdpdval(0))
        options%vdpddensityatvessel = 250.0_R8
        options%vdpddensityatxp     = 250.0_R8
        options%vdpddensityatinf    = 10.0_R8

        ! Grid generation approach
        options%verbosity           = 1
        options%ggmethod            = 'independent'
        options%cellconstructionmethod  = 'quads_triangles'
        options%legalcellstyle      = 'graph' 
        options%TMcellgriddingorder = 'sequential'
        options%readexistingrefdata = .false. 
        options%refdatafile         = './output/refdataTM.dat'

        ! Poloidal refinement options ('lengthbased' refinement options only)
        options%refmeth         = 'no'      
        options%reflengthtype   = 'euler'   
        options%refLBdoxp       = .true. 
        options%refLBdovessel   = .false. 
        options%solpstargetrefinements = .false. 
        options%refLBLmininf    = 0.0_R8
        options%refLBLmaxinf    = 100 ! some absurd big number
        options%refLBLminxp     = 0.0_R8
        options%refLBLmaxxp     = 100 ! some absurd big number  
        options%refLBdecaylengthxp = 0.1 
        allocate(options%refLBdecaylengthstructure(0), &
            options%refLBdecaylengthvert(0), options%refLBstructureIDs(0), &
            options%refLBvertIDs(0), options%refLBLminstructure(0), &
            options%refLBLmaxstructure(0))
        
        ! Radial refinement options
        options%radrefmeth         = 'no'   
        options%radreflengthtype   = 'euler'   
        options%radrefLBdosp       = .true. 
        options%radrefLBLmininf    = 0.0_R8
        options%radrefLBLmaxinf    = 100_R8 ! some absurd big number
        options%radrefLBLminsp     = 0.0_R8
        options%radrefLBLmaxsp     = 100_R8 ! some absurd big number  
        options%radrefLBdecaylengthsp = 0.1_R8 
        allocate(options%vdrdx(0), options%vdrdy(0), options%vdrdd(0), &
            options%vdrdval(0))

        ! Poloidal boundary layer options
        options%refBLdotarget   = .false. 
        options%refBLdovessel   = .false. 
        options%refBLdostructure    = .false. 
        options%refBLnctarget   = 0
        options%refBLncvessel   = 0
        options%refBLncstructure    = 0
        allocate(options%refBLdltarget(options%refBLnctarget), &
            options%refBLdlvessel(options%refBLncvessel), &
            options%refBLdlstructure(options%refBLncstructure), &
            options%refBLstructureID(0))
        options%refBLdltarget   = 0.001_R8 ! in m 
        options%refBLdlvessel   = 0.01_R8 ! in m
        options%refBLdlstructure    = 0.001_R8 ! in m

        ! Radial boundary layer options
        options%radrefBLdosp    = .false. 
        options%radrefBLncsp    = 0
        allocate(options%radrefBLdlsp(options%radrefBLncsp))
        options%radrefBLdlsp       = 0.001_R8 ! in m 

        ! Streamline tracer options
        options%orthtracerstep = 0.5
        options%orthtracernsteps = 2000

        ! Options for radial vertex distribution
        options%vdrtype             = 'uniform'
        options%vdrdfacelength      = 0.01 
        options%vdrdfieldwidth      = 4e-3
        options%vdrddecaylengthxp   = 0.005
        options%vdrddensityatxp     = 2500.0_R8
        options%vdrddensityatinf    = 250.0_R8
        options%vdrddecaylengthsp   = 0.1_R8
        options%vdrddensityatsp     = 1000.0_R8
        options%vdrdoxp             = .false.
        options%vdrdosp             = .true.

        ! Options for flux tube extensions
        options%extendtptubes       = .true. 
        options%extendvesseltubes   = .false. 
        options%evtmaxvessellength  = 0.2
        options%evtnoBL             = .true.

        ! Options for flux surface removal 
        options%removefluxsurfaces = .true.
        options%remfspsitol = 1e-4 
        options%remfspsirattol = 1e-1 
        
        ! Options for boundary triangle removal
        options%removenarrowboundarytriangles = .true. 
        options%rembndtriacriterion = 'angle' 
        options%rembndtriaminangle = 15.0_R8*pi_R8/180.0_R8

        ! Options for small face removal
        options%removefaces = .false. 
        options%remfacescriterion = 'facelength_radial_bnd'
        options%remfacesminlength = 2e-3        

        ! Diagnostics
        options%dogriddiagnostics = .true. ! default true 

        ! Label translation
        options%structurebasedlabels    = .false.
        options%forceSOLPStopology      = .false.
        options%SOLPStopology           = ''


    end subroutine 

    !------------------------------------------------------------------!
    !                            Option readers                        !
    !------------------------------------------------------------------!

    ! Goat options reader
    subroutine ReadGoatOptions(options)

        ! Description
        !============
        ! This routine reads in the goat options from a file of which 
        ! the full path should be given in options%inputfilepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown. 

        ! Declare variables
        !==================
        ! Arguments
        class(GoatoptionsUDT)            :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof, throwerror

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 
        throwerror = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGoatOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGoatOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        field = 'goat.debug'
        call ExtractOptionValueLogical0D(fid, field, options%debug)
        field = 'goat.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)

        ! Input filenames
        field = 'goat.gridfilepath'
        call ExtractOptionValueCharacter(fid, field, options%gridfilepath)
        field = 'goat.structurefilepath'
        call ExtractOptionValueCharacter(fid, field, options%structurefilepath)
        field = 'goat.magneticfieldfilepath'
        call ExtractOptionValueCharacter(fid, field, options%magneticfieldfilepath)
        field = 'goat.GDinputfilepath'
        call ExtractOptionValueCharacter(fid, field, options%gdinputfilepath)

        ! Output options
        field = 'goat.writefilepath'
        call ExtractOptionValueCharacter(fid, field, options%writefilepath)
        field = 'goat.write_final'
        call ExtractOptionValueLogical0D(fid, field, options%write_final)
        field = 'goat.write_traduitb2us'
        call ExtractOptionValueLogical0D(fid, field, options%write_traduitb2us)
        field = 'goat.write_b2agdat'
        call ExtractOptionValueLogical0D(fid, field, options%write_b2agdat)
        field = 'goat.write_OMPdata'
        call ExtractOptionValueLogical0D(fid, field, options%write_OMPdata)
        field = 'goat.write_topologicaldata'
        call ExtractOptionValueLogical0D(fid, field, options%write_topologicaldata)

        ! Case identification options
        field = 'goat.vesselmode'
        call ExtractOptionValueLogical0D(fid, field, options%vesselmode)
        field = 'goat.slab'
        call ExtractOptionValueLogical0D(fid, field, options%slab)
        field = 'goat.artificial_slab'
        call ExtractOptionValueLogical0D(fid, field, options%artificial_slab)

        ! Face label mappings
        field = 'goat.facelabelmappingGG'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelmappingGG)
        field = 'goat.facelabelmappingGD'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelmappingGD) 
        field = 'goat.facelabelmappingGA'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelmappingGA) 
        field = 'goat.facelabelsubfrom'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelsubfrom)
        field = 'goat.facelabelsubto'
        call ExtractOptionValueInteger1D(fid, field, &
            options%facelabelsubto)

        ! OMP and IMP
        field = 'goat.OMP_r'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_r)
        field = 'goat.OMP_z'
        call ExtractOptionValueReal1D(fid, field, &
            options%OMP_z)
        field = 'goat.IMP_r'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_r)
        field = 'goat.IMP_z'
        call ExtractOptionValueReal1D(fid, field, &
            options%IMP_z)
        
        ! Checks
        !=======
        ! Check for name clashes when using solps
        if (solps) then 
            ! Check input file names and assert they are equal to 
            ! assumed solps filenames
            if (options%gridfilepath /= solps_gridfilepath) then 
                print *, 'ReadGoatOptions: set goat.gridfilepath value ' // &
                    'to "' // solps_gridfilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%magneticfieldfilepath /= solps_magneticfieldfilepath) then 
                print *, 'ReadGoatOptions: set goat.magneticfieldfilepath value ' // &
                    'to "' // solps_magneticfieldfilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%structurefilepath /= solps_structurefilepath) then 
                print *, 'ReadGoatOptions: set goat.structurefilepath value ' // &
                    'to "' // solps_structurefilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 
            if (options%writefilepath /= solps_writefilepath) then 
                print *, 'ReadGoatOptions: set goat.writefilepath value ' // &
                    'to "' // solps_writefilepath // '" in GOAToptions.dat'
                throwerror = .true. 
            end if 

            ! Throw error
            if (throwerror) then 
                call gdErrorHandler('ReadGoatOptions: exiting due to ' // &
                    'wrong I/O filenames for SOLPS (see messages above)')
            end if

        end if 

        ! Other checks
        if (size(options%facelabelmappingGG)/=size(options%facelabelmappingGA)) then 
            call gdErrorHandler('ReadGOAToptions: facelabelmapping has inconsistent lengths')
        end if 

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Grid adaptation options reader
    subroutine ReadGAOptions(options)

        ! Description
        !============
        ! This routine reads in the ga options from a file of which 
        ! the full path should be given in options%inputfilepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown.
        
        
        ! Declare variables
        !==================
        ! Arguments
        class(GAoptionsUDT)            :: options  

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof, throwerror

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 
        throwerror = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGAOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGAOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if  
        
        ! Read options
        !=============
        ! General
        field = 'ga.plt'
        call ExtractOptionValueLogical0D(fid, field, options%plt)
        field = 'ga.plt_qm'
        call ExtractOptionValueLogical0D(fid, field, options%plt_qm)              
        field = 'ga.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth) 
        
        ! Operation options
        !==================

        ! Small triangles
        field = 'ga.rem_small_trias'
        call ExtractOptionValueLogical0D(fid, field, options%rem_small_trias)        
        field = 'ga.cut_off_pol'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_pol) 
        field = 'ga.cut_off_surf'
        call ExtractOptionValueReal0D(fid, field, options%cut_off_surf)

        ! Stacked triangles
        field = 'ga.stacked_trias'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_trias)          
        field = 'ga.stacked_trias_checkAR'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_trias_checkAR)     
        field = 'ga.stacked_trias_maxAR'
        call ExtractOptionValueReal0D(fid, field, options%stacked_trias_maxAR)
        field = 'ga.merge_stacked_trias'
        call ExtractOptionValueLogical0D(fid, field, options%merge_stacked_trias)
        field = 'ga.merge_stacked_trias_angle_threshold'
        call ExtractOptionValueReal0D(fid, field, options%merge_stacked_trias_angle_threshold)
        field = 'ga.merge_trap_into_stacked'
        call ExtractOptionValueLogical0D(fid, field, options%merge_trap_into_stacked) 
        
        ! Stacked to cutcell
        field = 'ga.stacked_to_cutcell'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_to_cutcell)
        field = 'ga.stacked_to_cutcell_nonuniform'
        call ExtractOptionValueLogical0D(fid, field, options%stacked_to_cutcell_nonuniform) 
               
        ! Splitting and merging
        field = 'ga.splitting'
        call ExtractOptionValueInteger1D(fid, field, options%splitting_array)    
        field = 'ga.pents_to_tria' 
        call ExtractOptionValueLogical0D(fid, field, options%pents_to_tria)
        field = 'ga.merging'
        call ExtractOptionValueInteger1D(fid, field, options%merging_array)
        field = 'ga.h_rad_threshold'
        call ExtractOptionValueReal0D(fid, field, options%h_rad_threshold)                
        field = 'ga.h_rad_core_threshold'
        call ExtractOptionValueReal0D(fid, field, options%h_rad_core_threshold)  

        ! Boundary layer grid
        field = 'ga.BLG'
        call ExtractOptionValueLogical0D(fid, field, options%BLG)
        field = 'ga.BLG_n_layers'
        call ExtractOptionValueInteger0D(fid, field, options%BLG_n_layers)
        field = 'ga.BLG_rescaling_factor'
        call ExtractOptionValueReal0D(fid, field, options%BLG_rescaling_factor)
        field = 'ga.BLG_smoothing_factors'
        call ExtractOptionValueReal1D(fid, field, options%BLG_smoothing_factors)

        ! Special operations
        field = 'ga.rem_stickout_trias'
        call ExtractOptionValueLogical0D(fid, field, options%rem_stickout_trias)        
        field = 'ga.rem_trias_tube'
        call ExtractOptionValueLogical0D(fid, field, options%rem_trias_tube)
        field = 'ga.rem_outershell'   
        call ExtractOptionValueLogical0D(fid, field, options%rem_outershell)
        field = 'ga.rem_tube_outershell_threshold'
        call ExtractOptionValueReal0D(fid, field, options%rem_tube_outershell_threshold)  
        field = 'ga.outershell_handling'
        call ExtractOptionValueCharacter(fid, field, options%outershell_handling)
        field = 'ga.rem_stickout_quad'        
        call ExtractOptionValueLogical0D(fid, field, options%rem_stickout_quad)
        field = 'ga.split_noalignedquads'
        call ExtractOptionValueLogical0D(fid, field, options%split_noalignedquads) 

        ! Splitting options
        !==================
        field = 'ga.no_pents'
        call ExtractOptionValueLogical0D(fid, field, options%no_pents)         
        field = 'ga.QTtype'
        call ExtractOptionValueCharacter(fid, field, options%QTtype) 
        field = 'ga.split_out'
        call ExtractOptionValueLogical0D(fid, field, options%split_out)
        field = 'ga.splittype'
        call ExtractOptionValueInteger1D(fid, field, options%splittype_array)
        field = 'ga.n_split'
        call ExtractOptionValueInteger1D(fid, field, options%n_split_array)
        field = 'ga.typeT'
        call ExtractOptionValueCharacter(fid, field, options%typeT)                                   
        field = 'ga.rad_type'
        call ExtractOptionValueInteger1D(fid, field, options%rad_type_array)                                   
        field = 'ga.pol_type'
        call ExtractOptionValueInteger1D(fid, field, options%pol_type_array)
        field = 'ga.dist_function_threshold_split'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_split)                                    
        field = 'ga.dist_function_threshold_split_wall'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_split_wall) 
        
        ! Merging options
        !================
        field = 'ga.no_hex'
        call ExtractOptionValueLogical0D(fid, field, options%no_hex)
        field = 'ga.merge_crit'
        call ExtractOptionValueInteger1D(fid, field, options%merge_crit_array) 
        field = 'ga.merge_h_pol_factor'
        call ExtractOptionValueReal0D(fid, field, options%merge_h_pol_factor)
        field = 'ga.n_merge'
        call ExtractOptionValueInteger1D(fid, field, options%n_merge_array)
        field = 'ga.merge_bias_limit'
        call ExtractOptionValueReal0D(fid, field, options%merge_bias_limit)        
        field = 'ga.dist_function_threshold_merge'
        call ExtractOptionValueReal0D(fid, field, options%dist_function_threshold_merge)  
        
        ! Pentagon options
        !=================
        field = 'ga.no_pents_area_merge'
        call ExtractOptionValueLogical0D(fid, field, options%no_pents_area_merge)        
        field = 'ga.no_pents_area_split'
        call ExtractOptionValueLogical0D(fid, field, options%no_pents_area_split)        
        field = 'ga.no_pents_area_type'
        call ExtractOptionValueCharacter(fid, field, options%no_pents_area_type)     
           
        field = 'ga.no_pents_area_maxR'
        call ExtractOptionValueReal0D(fid, field, options%no_pents_area_maxR)        
        field = 'ga.no_pents_area_minR'
        call ExtractOptionValueReal0D(fid, field, options%no_pents_area_minR)        
        field = 'ga.no_pents_area_maxZ'
        call ExtractOptionValueReal0D(fid, field, options%no_pents_area_maxZ)        
        field = 'ga.no_pents_area_minZ'
        call ExtractOptionValueReal0D(fid, field, options%no_pents_area_minZ)    
        
        ! Distance function options
        !==========================
        field = 'ga.dist_function'
        call ExtractOptionValueLogical0D(fid, field, options%dist_function)
        field = 'ga.base_func'
        call ExtractOptionValueCharacter(fid, field, options%base_func)        
        field = 'ga.d_rescale'     
        call ExtractOptionValueReal0D(fid, field, options%d_rescale)                             
        field = 'ga.d_rescale_wall'     
        call ExtractOptionValueReal0D(fid, field, options%d_rescale_wall)
        field = 'ga.dist_type'
        call ExtractOptionValueCharacter(fid, field, options%dist_type)
        field = 'ga.dist_type_wall'
        call ExtractOptionValueCharacter(fid, field, options%dist_type_wall)
        field = 'ga.d_char_type'
        call ExtractOptionValueCharacter(fid, field, options%d_char_type)
        field = 'ga.plt_dist_func'
        call ExtractOptionValueLogical0D(fid, field, options%plt_dist_func)

        ! Aposteriori
        !============
        field = 'ga.vxvol_style'
        call ExtractOptionValueInteger0D(fid, field, options%vxVol_style)
        field = 'ga.apost_interpolation_meth'
        call ExtractOptionValueCharacter(fid, field, options%apost_interpolation_meth)
        field = 'ga.apost_interpolationC'
        call ExtractOptionValueInteger0D(fid, field, options%apost_interpolationC)        
        field = 'ga.apost_interpolationM'
        call ExtractOptionValueInteger0D(fid, field, options%apost_interpolationM)                
        field = 'ga.apost_meth'
        call ExtractOptionValueCharacter(fid, field, options%apost_meth)
        field = 'ga.apost_use_na'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_na)
        field = 'ga.apost_use_ua'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_ua)
        field = 'ga.apost_use_te'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_te)
        field = 'ga.apost_use_ti'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_ti)
        field = 'ga.apost_use_tn'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_tn)
        field = 'ga.apost_use_po'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_po)
        field = 'ga.apost_use_kt'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_kt)
        field = 'ga.apost_use_zt'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_zt)
        field = 'ga.apost_use_resco'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_resco)
        field = 'ga.apost_use_resmo'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_resmo)
        field = 'ga.apost_use_resmt'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_resmt)
        field = 'ga.apost_use_reshe'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_reshe)
        field = 'ga.apost_use_reshi'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_reshi)
        field = 'ga.apost_use_reshn'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_reshn)
        field = 'ga.apost_use_respo'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_respo)
        field = 'ga.apost_use_reskt'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_reskt)
        field = 'ga.apost_use_reszt'
        call ExtractOptionValueLogical0D(fid, field, options%apost_use_reszt)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Grid deformation options reader
    subroutine ReadGDOptions(options)

        ! Description
        !============
        ! This routine reads in the grid deformation options from a file of which 
        ! the full path should be given in options%inputfilepath. The default
        ! options should have already been set at this point, as this 
        ! routine will only overwrite options that are present in the 
        ! user-specified options file. If no options file is present, 
        ! nothing is read in and a message will be shown. 

        ! Declare variables
        !==================
        ! Arguments
        class(GDoptionsUDT)             :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGDOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGDOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General
        field = 'gd.main.runtype'
        call ExtractOptionValueCharacter(fid, field, options%runtype)
        field = 'gd.main.gridtype'
        call ExtractOptionValueCharacter(fid, field, options%gridtype)
        field = 'gd.main.meth'
        call ExtractOptionValueCharacter(fid, field, options%meth)

        ! Files to read other options
        field = 'gd.main.designoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%designoptionsfile)
        field = 'gd.main.gridoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%gridoptionsfile)
        field = 'gd.main.magneticfieldoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%magneticfieldoptionsfile)
        field = 'gd.main.numparamsoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%numparamsoptionsfile)
        field = 'gd.main.environmentoptionsfile'
        call ExtractOptionValueCharacter(fid, field, options%environmentoptionsfile)
        
        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)


    end subroutine

    ! Grid options reader
    subroutine ReadGridOptions(options)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(GridOptionsUDT)            :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false.

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGridOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGridOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! Gridtype
        field = 'goat.grid.type'
        call ExtractOptionValueCharacter(fid, field, options%type)

        ! Reading method
        field = 'goat.grid.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Magnetic field options reader
    subroutine ReadMagneticFieldOptions(options)

        ! Description
        !============
        ! Read in grid options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(MagneticFieldOptionsUDT)            :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadMagneticFieldOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadMagneticFieldOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'goat.mf.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)

        ! Interpolant
        field = 'goat.mf.interpC'
        call ExtractOptionValueInteger0D(fid, field, options%interpC)
        field = 'goat.mf.interpM'
        call ExtractOptionValueInteger0D(fid, field, options%interpM)
        field = 'goat.mf.interpmeth'
        call ExtractOptionValueCharacter(fid, field, options%interpmeth)
        field = 'goat.mf.reinterpolate'
        call ExtractOptionValueLogical0D(fid, field, options%reinterpolate)
        field = 'goat.mf.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'goat.mf.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)

        ! RBtor
        field  = 'goat.mf.RBtor'
        call ExtractOptionValueReal0D(fid, field, options%RBtor)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Vessel options reader
    subroutine ReadVesselOptions(options)

        ! Description
        !============
        ! Read in vessel options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(VesselOptionsUDT)         :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadVesselOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadVesselOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'goat.vessel.readmeth'
        call ExtractOptionValueCharacter(fid, field, options%readmeth)
        field = 'goat.vessel.exclude'
        call ExtractOptionValueInteger1D(fid, field, options%exclude)

        ! Refinement
        field = 'goat.vessel.refinevessel'
        call ExtractOptionValueLogical0D(fid, field, options%refine)
        field = 'goat.vessel.maxvesseldist'
        call ExtractOptionValueReal0D(fid, field, options%maxdist)
        field = 'goat.vessel.minreffac'
        call ExtractOptionValueInteger0D(fid, field, options%minreffac)
        
        ! Target plates
        field = 'goat.vessel.TP'
        call ExtractOptionValueInteger1D(fid, field, options%TP)
        field = 'goat.vessel.TPind'
        call ExtractOptionValueInteger1D(fid, field, options%TPind)
        field = 'goat.vessel.solps_limiter_target'
        if (.not. allocated(options%solpslimitertarget)) then
            allocate(options%solpslimitertarget(0))
        end if
        call ExtractOptionValueInteger1D(fid, field, &
            options%solpslimitertarget)
        field = 'goat.vessel.solps_lower_divertor_inner_target'
        if (.not. allocated(options%solpslowdivinnertarget)) then
            allocate(options%solpslowdivinnertarget(0))
        end if
        call ExtractOptionValueInteger1D(fid, field, &
            options%solpslowdivinnertarget)
        field = 'goat.vessel.solps_lower_divertor_outer_target'
        if (.not. allocated(options%solpslowdivoutertarget)) then
            allocate(options%solpslowdivoutertarget(0))
        end if
        call ExtractOptionValueInteger1D(fid, field, &
            options%solpslowdivoutertarget)
        field = 'goat.vessel.solps_upper_divertor_inner_target'
        if (.not. allocated(options%solpsupdivinnertarget)) then
            allocate(options%solpsupdivinnertarget(0))
        end if
        call ExtractOptionValueInteger1D(fid, field, &
            options%solpsupdivinnertarget)
        field = 'goat.vessel.solps_upper_divertor_outer_target'
        if (.not. allocated(options%solpsupdivoutertarget)) then
            allocate(options%solpsupdivoutertarget(0))
        end if
        call ExtractOptionValueInteger1D(fid, field, &
            options%solpsupdivoutertarget)


        ! Shape representation options
        field = 'goat.vessel.shapemeth'
        call ExtractOptionValueCharacter(fid, field, options%shapemeth)
        field = 'goat.vessel.resx'
        call ExtractOptionValueInteger0D(fid, field, options%resx)
        field = 'goat.vessel.resy'
        call ExtractOptionValueInteger0D(fid, field, options%resy)
        field = 'goat.vessel.interpC'
        call ExtractOptionValueInteger0D(fid, field, options%interpC)
        field = 'goat.vessel.interpM'
        call ExtractOptionValueInteger0D(fid, field, options%interpM)
        field = 'goat.vessel.offsetfracx'
        call ExtractOptionValueReal0D(fid, field, options%offsetfracx)
        field = 'goat.vessel.offsetfracy'
        call ExtractOptionValueReal0D(fid, field, options%offsetfracy)
        field = 'goat.vessel.xrange'
        call ExtractOptionValueReal1D(fid, field, options%xrange)
        field = 'goat.vessel.yrange'
        call ExtractOptionValueReal1D(fid, field, options%yrange)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Environment options reader
    subroutine ReadEnvironmentOptions(options)

        ! Description
        !============
        ! Read in environment options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(EnvironmentOptionsUDT)    :: options 

        ! Auxiliary
        integer                         :: openstatus
        character(:), allocatable       :: field 
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadEnvironmentOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadEnvironmentOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        field = 'goat.statefilepath'
        call ExtractOptionValueCharacter(fid, field, options%statefilepath)        
        field = 'goat.readstate'
        call ExtractOptionValueLogical0D(fid, field, options%readstate)
        field = 'goat.readstatemeth'
        call ExtractOptionValueCharacter(fid, field, options%readstatemeth)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Topomesh options reader
    subroutine ReadTopomeshOptions(options)

        ! Description
        !============
        ! Read in Topomesh options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(TopomeshOptionsUDT)       :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadTopomeshOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadTopomeshOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! I/O
        field = 'gg.tm.readexistingTM'
        call ExtractOptionValueLogical0D(fid, field, options%readexistingTM)
        field = 'gg.tm.readexistingtracers'
        call ExtractOptionValueLogical0D(fid, field, options%readexistingtracers)
        field = 'gg.tm.writedebugoutput'
        call ExtractOptionValueLogical0D(fid, field, options%writedebugoutput)
        field = 'gg.tm.TMfilepath'
        call ExtractOptionValueCharacter(fid, field, options%TMfilepath)
        field = 'gg.tm.TMfieldtracerfilepath'
        call ExtractOptionValueCharacter(fid, field, options%TMfieldtracerfilepath)
        field = 'gg.tm.TMvesseltracerfilepath'
        call ExtractOptionValueCharacter(fid, field, options%TMvesseltracerfilepath)

        ! Resolution 
        field = 'gg.tm.field.resx'
        call ExtractOptionValueInteger0D(fid, field, options%fresx)
        field = 'gg.tm.field.resy'
        call ExtractOptionValueInteger0D(fid, field, options%fresy)
        field = 'gg.tm.vessel.resx'
        call ExtractOptionValueInteger0D(fid, field, options%vresx)
        field = 'gg.tm.vessel.resy'
        call ExtractOptionValueInteger0D(fid, field, options%vresy)
        field = 'gg.tm.contour.npmin'
        call ExtractOptionValueInteger0D(fid, field, options%npmin)
        field = 'gg.tm.contour.npmax'
        call ExtractOptionValueInteger0D(fid, field, options%npmax)
        field = 'gg.tm.contour.dl'
        call ExtractOptionValueReal0D(fid, field, options%dl)

        ! Refinement
        field = 'gg.tm.field.donewton'
        call ExtractOptionValueLogical0D(fid, field, options%fdonewton)
        field = 'gg.tm.vessel.donewton'
        call ExtractOptionValueLogical0D(fid, field, options%vdonewton)
        field = 'gg.tm.fieldtol'
        call ExtractOptionValueReal0D(fid, field, options%ffieldtol)

        ! Additional options
        field = 'gg.tm.dotpvesselbased'
        call ExtractOptionValueLogical0D(fid, field, options%dotpvesselbased)

        ! Adaptations
        field = 'gg.tm.doadaptations'
        call ExtractOptionValueLogical0D(fid, field, options%doadaptations)

        field = 'gg.tm.addcoreboundaries'
        call ExtractOptionValueLogical0D(fid, field, options%addcoreboundaries)
        field = 'gg.tm.addPFboundaries'
        call ExtractOptionValueLogical0D(fid, field, options%addPFboundaries)
        field = 'gg.tm.removecoreregions'
        call ExtractOptionValueLogical0D(fid, field, options%removecoreregions)
        field = 'gg.tm.removewidegridregions'
        call ExtractOptionValueLogical0D(fid, field, options%removewidegridregions)
        field = 'gg.tm.cbnd.frac'
        call ExtractOptionValueReal0D(fid, field, options%coreboundariesfrac)
        field = 'gg.tm.PFbnd.frac'
        call ExtractOptionValueReal0D(fid, field, options%PFboundariesfrac)
        field = 'gg.tm.removenoncoreregions'
        call ExtractOptionValueLogical0D(fid, field, options%removenoncoreregions)
        field = 'gg.tm.mergetubes'
        call ExtractOptionValueLogical0D(fid, field, options%mergetubes)
        field = 'gg.tm.mergetubes.meth'
        call ExtractOptionValueCharacter(fid, field, options%mergetubemeth)
        field = 'gg.tm.mergetubes.allowseparatrix'
        call ExtractOptionValueLogical0D(fid, field, options%mtallowseparatrix)
        field = 'gg.tm.mergetubes.allowcore'
        call ExtractOptionValueLogical0D(fid, field, options%mtallowcore)
        field = 'gg.tm.mergetubes.allowpf'
        call ExtractOptionValueLogical0D(fid, field, options%mtallowpf)
        field = 'gg.tm.dpsimintubes'
        call ExtractOptionValueReal0D(fid, field, options%dpsimintubes)
        field = 'gg.tm.lradmintubes'
        call ExtractOptionValueReal0D(fid, field, options%lradmintubes)
        field = 'gg.tm.mergetubes.oldstyle'
        call ExtractOptionValueLogical0D(fid, field, options%mtoldstyle)
        field = 'gg.tm.lradmintubes'
        call ExtractOptionValueReal0D(fid, field, options%lradmintubes)
        field = 'gg.tm.simplifyseparatrix'
        call ExtractOptionValueLogical0D(fid, field, options%simplifyseparatrix)
        field = 'gg.tm.simplifyseparatrix.lpolsepmin'
        call ExtractOptionValueReal0D(fid, field, options%lpolsepmin)

        field = 'gg.tm.removevesselregions'
        call ExtractOptionValueLogical0D(fid, field, options%removevesselregions)
        field = 'gg.tm.rvrvesselIDs'
        call ExtractOptionValueInteger1D(fid, field, options%rvrvesselIDs)
        field = 'gg.tm.rvrfullycovered'
        call ExtractOptionValueLogical0D(fid, field, options%rvrfullycovered)
        field = 'gg.tm.rvrretain'
        call ExtractOptionValueLogical0D(fid, field, options%rvrretain)
        field = 'gg.tm.rvrdocascade'
        call ExtractOptionValueLogical0D(fid, field, options%rvrdocascade)
        field = 'gg.tm.rvrcascadedir'
        call ExtractOptionValueCharacter(fid, field, options%rvrcascadedir)

        field = 'gg.tm.alignvesselparts'
        call ExtractOptionValueLogical0D(fid, field, options%alignvesselparts)
        field = 'gg.tm.avprefinevessel'
        call ExtractOptionValueLogical0D(fid, field, options%avprefinevessel)
        field = 'gg.tm.avpminangle'
        call ExtractOptionValueReal0D(fid, field, options%avpminangle)
        field = 'gg.tm.avpmaxvesseldist'
        call ExtractOptionValueReal0D(fid, field, options%avpmaxvesseldist)
        field = 'gg.tm.avpminreffac'
        call ExtractOptionValueInteger0D(fid, field, options%avpminreffac)


        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

    ! Grid generator options reader
    subroutine ReadGGOptions(options)

        ! Description
        !============
        ! Read in GG options from file. It is assumed that the 
        ! filepath has been set correctly. 

        ! Declare variables
        !==================
        ! Arguments
        class(GGOptionsUDT)             :: options 

        ! Auxiliary
        integer                         :: openstatus 
        character(:), allocatable       :: field
        integer, parameter              :: fid = 10 
        logical                         :: reachedeof

        ! Initialize
        !===========
        ! Variables
        reachedeof = .false. 

        ! Open the file, check if it exists
        open(unit=fid, file=options%inputfilepath, status='old', &
            iostat=openstatus)

        if (openstatus > 0) then 
            ! Something wrong when reading file - continue with default
            ! values
            print *, 'ReadGGOptions: could not open file, ' &
                // 'taking default options...'
        elseif (openstatus < 0) then 
            ! File appears to be empty
            print *, 'ReadGGOptions: file appears to be empty, ' &
                // 'taking default options...'
        end if
        
        ! Read options
        !=============
        ! General options
        field  = 'gg.verbosity'
        call ExtractOptionValueInteger0D(fid, field, options%verbosity)
        field  = 'gg.vertexconstructionmethod'
        call ExtractOptionValueCharacter(fid, field, options%ggmethod)
        field  = 'gg.cellconstructionmethod'
        call ExtractOptionValueCharacter(fid, field, options%cellconstructionmethod)
        field  = 'gg.legalcellstyle'
        call ExtractOptionValueCharacter(fid, field, options%legalcellstyle)
        field = 'gg.TMcellgriddingorder'
        call ExtractOptionValueCharacter(fid, field, options%TMcellgriddingorder)

        ! Label translation
        field = 'gg.labels.structurebased'
        call ExtractOptionValueLogical0D(fid, field, options%structurebasedlabels)
        field = 'gg.labels.forceSOLPStopology'
        call ExtractOptionValueLogical0D(fid, field, options%forceSOLPStopology)
        field = 'gg.labels.SOLPStopology'
        call ExtractOptionValueCharacter(fid, field, options%SOLPStopology)

        ! Refinement options (general)
        field = 'gg.ref.meth'
        call ExtractOptionValueCharacter(fid, field, options%refmeth) 
        field = 'gg.ref.LB.lengthtype'
        call ExtractOptionValueCharacter(fid, field, options%reflengthtype)
        field = 'gg.ref.readexistingrefdata'
        call ExtractOptionValueLogical0D(fid, field, options%readexistingrefdata) 
        field = 'gg.ref.refdatafile'
        call ExtractOptionValueCharacter(fid, field, options%refdatafile) 

        ! Length-based refinement options (poloidal)
        field  = 'gg.ref.LB.doxp'
        call ExtractOptionValueLogical0D(fid, field, options%refLBdoxp)
        field  = 'gg.ref.LB.dovessel'
        call ExtractOptionValueLogical0D(fid, field, options%refLBdovessel)
        field  = 'gg.ref.LB.Lmininf'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmininf)
        field  = 'gg.ref.LB.Lmaxinf'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmaxinf)
        
        field  = 'gg.ref.LB.Lminxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBLminxp)
        field  = 'gg.ref.LB.Lmaxxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBLmaxxp)
        field  = 'gg.ref.LB.decaylengthxp'
        call ExtractOptionValueReal0D(fid, field, options%refLBdecaylengthxp)

        field  = 'gg.ref.LB.Lminstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBLminstructure)
        field  = 'gg.ref.LB.Lmaxstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBLmaxstructure)
        field  = 'gg.ref.LB.decaylengthstructure'
        call ExtractOptionValueReal1D(fid, field, options%refLBdecaylengthstructure)
        field  = 'gg.ref.solps_target_refinements'
        call ExtractOptionValueLogical0D(fid, field, &
            options%solpstargetrefinements)
        field  = 'gg.ref.LB.structureIDs'   
        call ExtractOptionValueInteger1D(fid, field, options%refLBstructureIDs)
        field  = 'gg.ref.LB.vertIDs'   
        call ExtractOptionValueInteger1D(fid, field, options%refLBvertIDs)

        ! Refinement options (radial)
        field = 'gg.radref.meth'
        call ExtractOptionValueCharacter(fid, field, options%radrefmeth) 
        field  = 'gg.radref.LB.lengthtype'
        call ExtractOptionValueCharacter(fid, field, options%radreflengthtype)
        field  = 'gg.radref.LB.dosp'
        call ExtractOptionValueLogical0D(fid, field, options%radrefLBdosp)
        field  = 'gg.radref.LB.Lmininf'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmininf)
        field  = 'gg.radref.LB.Lmaxinf'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmaxinf)
        
        field  = 'gg.radref.LB.Lminsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLminsp)
        field  = 'gg.radref.LB.Lmaxsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBLmaxsp)
        field  = 'gg.radref.LB.decaylengthsp'
        call ExtractOptionValueReal0D(fid, field, options%radrefLBdecaylengthsp)

        ! Boundary layer options (only for length-based ref, poloidal)
        field = 'gg.ref.BL.dotarget'
        call ExtractOptionValueLogical0D(fid, field, options%refBLdotarget)
        field = 'gg.ref.BL.dovessel'
        call ExtractOptionValueLogical0D(fid, field, options%refBLdovessel)
        field = 'gg.ref.BL.dostructure'
        call ExtractOptionValueLogical0D(fid, field, options%refBLdostructure)
        field = 'gg.ref.BL.nctarget'
        call ExtractOptionValueInteger0D(fid, field, options%refBLnctarget)
        field = 'gg.ref.BL.ncvessel'
        call ExtractOptionValueInteger0D(fid, field, options%refBLncvessel)
        field = 'gg.ref.BL.ncstructure'
        call ExtractOptionValueInteger0D(fid, field, options%refBLncstructure)
        field = 'gg.ref.BL.dltarget'
        call ExtractOptionValueReal1D(fid, field, options%refBLdltarget)
        field = 'gg.ref.BL.dlvessel'
        call ExtractOptionValueReal1D(fid, field, options%refBLdlvessel)
        field = 'gg.ref.BL.dlstructure'
        call ExtractOptionValueReal1D(fid, field, options%refBLdlstructure)
        field = 'gg.ref.BL.dllengthbased'
        call ExtractOptionValueLogical0D(fid, field, options%refdlBLlengthbased)
        field = 'gg.ref.BL.structureID'
        call ExtractOptionValueInteger1D(fid, field, options%refBLstructureID)

        ! Boundary layer options (only for length-based ref, radial)
        field = 'gg.radref.BL.dosp'
        call ExtractOptionValueLogical0D(fid, field, options%radrefBLdosp)
        field = 'gg.radref.BL.ncsp'
        call ExtractOptionValueInteger0D(fid, field, options%radrefBLncsp)
        field = 'gg.radref.BL.dlsp'
        call ExtractOptionValueReal1D(fid, field, options%radrefBLdlsp)
        field = 'gg.radref.BL.dllengthbased'
        call ExtractOptionValueLogical0D(fid, field, options%radrefdlBLlengthbased)

        ! Contouring options in grid generator
        field = 'gg.vd.contouring.resx'
        call ExtractOptionValueInteger0D(fid, field, options%gcresx)
        field = 'gg.vd.contouring.resy'
        call ExtractOptionValueInteger0D(fid, field, options%gcresy)
        field = 'gg.vd.contouring.coarsen'
        call ExtractOptionValueLogical0D(fid, field, options%coarsencontours)

        ! Orthogonal line tracer options
        field = 'gg.vd.orthlinetracing.step'
        call ExtractOptionValueReal0D(fid, field, options%orthtracerstep)
        field = 'gg.vd.orthlinetracing.nsteps'
        call ExtractOptionValueInteger0D(fid, field, options%orthtracernsteps)

        ! Options for poloidal vertex distribution
        field = 'gg.vd.pd.type'
        call ExtractOptionValueCharacter(fid, field, options%vdptype)
        field = 'gg.vd.pd.distribution.plftype'
        call ExtractOptionValueCharacter(fid, field, options%vdpplftype)
        field = 'gg.vd.pd.distribution.facelength'
        call ExtractOptionValueReal0D(fid, field, options%vdpdfacelength)
        field = 'gg.vd.pd.distribution.decaylengthplf'
        call ExtractOptionValueReal0D(fid, field, options%vdpddecaylengthplf)
        field = 'gg.vd.pd.distribution.decaylengthxp'
        call ExtractOptionValueReal0D(fid, field, options%vdpddecaylengthxp)
        field = 'gg.vd.pd.distribution.vdpdincludexp'
        call ExtractOptionValueLogical0D(fid, field, options%vdpdincludexp)
        field = 'gg.vd.pd.distribution.points.x'
        call ExtractOptionValueReal1D(fid, field, options%vdpdx)
        field = 'gg.vd.pd.distribution.points.y'
        call ExtractOptionValueReal1D(fid, field, options%vdpdy)
        field = 'gg.vd.pd.distribution.points.d'
        call ExtractOptionValueReal1D(fid, field, options%vdpdd)
        field = 'gg.vd.pd.distribution.points.val'
        call ExtractOptionValueReal1D(fid, field, options%vdpdval)
        field = 'gg.vd.pd.distribution.densityatvessel'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatvessel)
        field = 'gg.vd.pd.distribution.densityatxp'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatxp)
        field = 'gg.vd.pd.distribution.densityatinf'
        call ExtractOptionValueReal0D(fid, field, options%vdpddensityatinf)

        ! Options for radial vertex distribution
        field = 'gg.vd.rd.type'
        call ExtractOptionValueCharacter(fid, field, options%vdrtype)
        field = 'gg.vd.rd.distribution.doxp'
        call ExtractOptionValueLogical0D(fid, field, options%vdrdoxp)
        field = 'gg.vd.rd.distribution.dosp'
        call ExtractOptionValueLogical0D(fid, field, options%vdrdosp)
        field = 'gg.vd.rd.distribution.fieldwidth'
        call ExtractOptionValueReal0D(fid, field, options%vdrdfieldwidth)
        field = 'gg.vd.rd.distribution.facelength'
        call ExtractOptionValueReal0D(fid, field, options%vdrdfacelength)
        field = 'gg.vd.rd.distribution.decaylengthxp'
        call ExtractOptionValueReal0D(fid, field, options%vdrddecaylengthxp)
        field = 'gg.vd.rd.distribution.densityatxp'
        call ExtractOptionValueReal0D(fid, field, options%vdrddensityatxp)
        field = 'gg.vd.rd.distribution.decaylengthsp'
        call ExtractOptionValueReal0D(fid, field, options%vdrddecaylengthsp)
        field = 'gg.vd.rd.distribution.densityatsp'
        call ExtractOptionValueReal0D(fid, field, options%vdrddensityatsp)
        field = 'gg.vd.rd.distribution.densityatinf'
        call ExtractOptionValueReal0D(fid, field, options%vdrddensityatinf)
        field = 'gg.vd.rd.distribution.points.x'
        call ExtractOptionValueReal1D(fid, field, options%vdrdx)
        field = 'gg.vd.rd.distribution.points.y'
        call ExtractOptionValueReal1D(fid, field, options%vdrdy)
        field = 'gg.vd.rd.distribution.points.d'
        call ExtractOptionValueReal1D(fid, field, options%vdrdd)
        field = 'gg.vd.rd.distribution.points.val'
        call ExtractOptionValueReal1D(fid, field, options%vdrdval)

        ! Options for extending flux tubes
        field = 'gg.adap.extendtptubes'
        call ExtractOptionValueLogical0D(fid, field, options%extendtptubes)
        field = 'gg.adap.extendvesseltubes'
        call ExtractOptionValueLogical0D(fid, field, options%extendvesseltubes)
        field = 'gg.adap.evt.maxvessellength'
        call ExtractOptionValueReal0D(fid, field, options%evtmaxvessellength)
        field = 'gg.adap.evt.noBL'
        call ExtractOptionValuelogical0D(fid, field, options%evtnoBL)

        ! Options for flux surface removal 
        field = 'gg.vd.removefluxsurfaces'
        call ExtractOptionValueLogical0D(fid, field, options%removefluxsurfaces)
        field = 'gg.vd.rfs.mark.psitol'
        call ExtractOptionValueReal0D(fid, field, options%remfspsitol)
        field = 'gg.vd.rfs.mark.psirattol'
        call ExtractOptionValueReal0D(fid, field, options%remfspsirattol)
        
        ! Options for boundary triangle removal
        field = 'gg.vd.removenarrowboundarytriangles'
        call ExtractOptionValueLogical0D(fid, field, options%removenarrowboundarytriangles)
        field = 'gg.vd.rnbt.mark.criterion'
        call ExtractOptionValueCharacter(fid, field, options%rembndtriacriterion)
        field = 'gg.vd.rnbt.mark.minangle'
        call ExtractOptionValueReal0D(fid, field, options%rembndtriaminangle)

        ! Options for small face removal
        field = 'gg.vd.removefaces'
        call ExtractOptionValueLogical0D(fid, field, options%removefaces)
        field = 'gg.vd.rf.criterion'
        call ExtractOptionValueCharacter(fid, field, options%remfacescriterion)
        field = 'gg.vd.rf.minlength'
        call ExtractOptionValueReal0D(fid, field, options%remfacesminlength)  
        
        ! Diagnostic options
        field = 'gg.dogriddiagnostics'
        call ExtractOptionValueLogical0D(fid, field, options%dogriddiagnostics)

        ! Housekeeping
        !=============
        ! Close the file
        close(unit=fid)

    end subroutine

end module goatmod_userinput
