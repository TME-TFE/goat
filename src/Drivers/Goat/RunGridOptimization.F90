subroutine RunGridOptimization(grid, magneticField, environment, &
    options)
    ! Description
    !============
    ! Driver for optimization-based grid deformation. It is assumed that
    ! the initial given grid, magnetic field, and environment structure 
    ! are properly initialized and set up. This routine additionally
    ! sets up the optimization problem, for which options should be 
    ! specified in the designoptions structure. Filepaths for the 
    ! designoptions user input file should be present in 
    ! options%designfilepath.

    ! Initialize
    !===========
    ! Declare modules
    use gdmod_types
    use gdmod_userinput 
    use gdmod_plots
    use mod_plotter
    use gdmod_optimizationengine

    ! The usual
    implicit none

    ! Declare variables
    type(GridUDT)                       :: grid
    type(MagneticFieldUDT)              :: magneticField
    type(OptimizationEngineGDUDT)       :: optimizationdriver
    type(GDoptionsUDT), intent(in)      :: options 
    type(DesignOptionsUDT)              :: designoptions
    ! type(NumOptionsUDT)             :: num
    type(EnvironmentUDT)                :: environment

    ! Debug
    logical                             :: makegridplots  = .false.

    character(:), allocatable           :: vertpath, cellpath

    ! Set additional options
    !=======================
    ! Set paths from where options should be read
    designoptions%inputfilepath = options%designoptionsfile

    ! Set optimization options
    call designoptions%Set()

    ! Set other numerical parameters
    !call SetNumOptions(num);

    ! Initialize
    !===========
    ! Initialize the grid design problem (as an optimization problem)
    call ConstructGridDesignProblem(optimizationdriver, designoptions, &
        grid, magneticField, environment)

    ! Write initial coordinates
    allocate(character(len('vertices_init')) :: vertpath)
    allocate(character(len('cells_init')) :: cellpath)
    vertpath = 'vertices_init'
    cellpath = 'cells_init'
    call WriteGridVertices(grid, vertpath)
    call WriteGridCells(grid, cellpath)

    ! Solve
    !======
    ! Simply call the solver ... 
    call optimizationdriver%Driver()

    ! Plot
    associate(problem => optimizationdriver%problem)
        
    select type (problem)

    type is (OptimizationProblemGDUDT)

        if (makegridplots) then 
            call PlotGridCells(problem%grid, '-p')
        end if 

        ! Unpack solution
        grid = problem%grid
        magneticField = problem%magneticField 
        environment = problem%environment

    end select

    end associate



end subroutine