# GOAT

Grid Generation, Optimization and Adapation Toolbox (acronyms: GOAT/Goat/goat).

## Description
Toolbox to generate (2D) magnetic field aligned grids for plasma edge simulations. Currently focused on the SOLPS plasma edge simulator, but in principle any solver can be coupled. Arbitrary magnetic field topologies and (closed) geometry are supported. Please cite the following papers when using:
- Grid generation: S. Van den Kerkhof, S. Carli, W. Dekeyser. "Automated 2D plasma edge grid generation for arbitrary magnetic topologies". Submitted for publication in Contributions to Plasma Physics
- Grid optimization/deformation: S. Van den Kerkhof, N. Vervloesem, S. Carli, W. Dekeyser. "Application of an automated grid deformation tool for divertor shape optimization in SOLPS-ITER", Contributions to Plasma Physics, 2024, 64(7-8), e202300134. https://doi.org/10.1002/ctpp.202300134
- Grid adaptation: N. Vervloesem, W. Dekeyser, S. Van den Kerkhof, M. Baelmans. "Error-based grid adaptation methods for plasma edge simulations with SOLPS-ITER", Contributions to Plasma Physics, 2024, 64(7-8), e202300126. https://doi.org/10.1002/ctpp.202300126

## Installation
To install goat, you need a Unix-like environment and the tcsh shell. Additionally, you need to have at least the following libraries and packages installed:
- Make 
- A fortran compiler that can cope with the F2008/F2018 standard (gfortran version 13.3.0 is currently being used by the developers)
- A C compiler (gcc version 13.3.0 is currently being used by the developers)
- SuiteSparse (see https://github.com/DrTimothyAldenDavis/SuiteSparse)
- OpenBLAS

Using Ubuntu, one can use the package manager to install these libraries (e.g. `sudo apt install libopenblas-dev` for OpenBLAS, `sudo apt install libsuitesparse-dev`, `sudo apt install gfortran`)
To run the Python visualization scripts, a Python installation with additional packages (and their dependencies) is required:
- Matplotlib
- Numpy
- Shapely

Having these dependencies installed, the code is installed as follows:
- open a clean tcsh terminal
- clone the git repository
- `cd` into the top folder after cloning and execute `source setup.csh`. This will load in the necessary environment variables for compilation.

### Using the old makefile (will become deprecated in the future)
- To compile goat, execute `make goat` (or `make goat_debug`) to compile the executable name `goat.exe` and `goat_debug.exe`. After compilation, these should be found in the executables again.
- To link properly to the newly compiled executables, rerun `source setup.csh` again

### Using CMake
- Running the 'build.sh' script in the main folder will build the basic 'Debug' and 'Release' versions of the code automatically. The resulting executables are stored in a different location and linked to using symbolic links 'goat_debug.exe' and 'goat.exe', such that, usage-wise, nothing should change (these paths are added in the setup.csh script, so executing 'goat_debug.exe' in a case folder should work exactly the same as before)
- Note that currently, any environment-specific compiler settings etc are ignored (unless specified in CMake variables before). This remains to be improved in the future. 

To test if installation was successful, try one of the example cases (see Usage section below). 

## Usage
To use goat, first compile it (see section Installation) to retrieve the goat executable `goat.exe` or `goat_debug.exe`. To run an example case (e.g. the ASDEX grid generation case), do the following steps:
- open a clean tcsh terminal, `cd` to the top goat directory (note: if you compiled for the first time and work in the same terminal, the shell cache may not have been updated with the newly constructed symlinks/executables. For tcsh, type `rehash` to refresh this, or open a new terminal)
- create a `setup.csh` file for your machine/cluster, if not already present
- execute `source setup.csh`
- go to the example case (i.e. execute `cd ./Examples/ASDEX`)
- run the example case by executing `goat.exe` or `goat_debug.exe`
- After succesful execution, a `traduit.out.b2us` file is produced, along with other data in `./output`. The grid is stored in the `traduit.out.b2us` file. For SOLPS, this is the input for the `b2ag` preprocessor (along with a `b2ag.dat` file)
- To visualize results, one can execute `pggoutput` in the terminal, which will read in the grid and additional data in the `output` folder (this may fail if non-standard names are used or if Python packages were not installed correctly)

To create a new case, create a new folder called `Runs` on the same level as `src`. This folder will be ignored by the current .gitignore settings. To run goat, only three inputs are required:
- a `GOAToptions.dat` file containing the user-specified options (best to copy-paste from an existing example)
- an `rzpsi.dat` or `.equ` file containing the magnetic equilibrium data in 2D
- a `structure.dat` file containing the vessel structures of the device. These structures should form a set of closed, nested, non-selfintersecting polygons, such that the interior of the domain is well defined.
To run this case, follow the steps above for the example case, but go `./Runs/<name_of_your_case>` when executing goat.

## Known issues
- The current make-based build system does not always capture all dependencies - improvements welcome! In case of issues, run `make deepclean` to remove all intermediate files and build folders such that recompilation happens from scratch.
- The Cmake approach currently only covers the goat executable, but not yet others (which are not typically used).

## Known behavior that is not an issue
- If one uses SOLPS and has sourced that setup file, the compilation and running of goat may be hampered. Solution: use a clean terminal as described above. For CMake, this shouldn't be a problem yet, since executables that would use SOLPS are not yet supported and implemented. 

## Support
In case of issues, please contact either the developer(s) at KU Leuven or their supervisors:
- Sander Van den Kerkhof (sander.vandenkerkhof@kuleuven.be) -> main developer
- Wouter Dekeyser (wouter.dekeyser@kuleuven.be) -> supervisor

## Contributing
Contributions by the community, be it bugfixes, build improvements, or completely new features, are very much welcome! For large or fundamental changes, please first contact the developers to discuss this further. Other improvements can be done through pull requests after review by the main devs. For more detailed guidelines on how to use the build system etc, please see the Contributing.md document. 

## Authors and acknowledgment
Please see the description for authors and reference papers for the different modules of goat. The authors are very grateful for the many interesting testcases provided by the plasma edge modelling community that helped shape this toolbox. 

## License
Goat is licensed under a EUPLv1.2 license (see also LICENSE.md). 

## Project status
Currently being actively developed and supported at KU Leuven. 
