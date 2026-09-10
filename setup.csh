#! /bin/tcsh -f

# Preliminary setup script for running goat in standalone mode (for 
# running as a submodule, see the overarching module setup workflow)
# Very rudimentary and simple, but should get the job done. 
# The following actions are done:
# - the path to 'goatrun' is added to the $PATH environment variable. 
# 'goatrun' is the main running program.
# - several environment variables are set to default values

# Say hi to the user
echo ' Welcome to the ... '
echo '   ______     _     __   ______                           __  _           '
echo '  / ____/____(_)___/ /  / ____/__  ____  ___  _________ _/ /_(_)___  ____ '
echo ' / / __/ ___/ / __  /  / / __/ _ \/ __ \/ _ \/ ___/ __ `/ __/ / __ \/ __ \'
echo '/ /_/ / /  / / /_/ /  / /_/ /  __/ / / /  __/ /  / /_/ / /_/ / /_/ / / / /'
echo '\____/_/  /_/\__,_/   \____/\___/_/ /_/\___/_/   \__,_/\__/_/\____/_/ /_/ '
echo '   ____        __  _           _             __  _                ___     '
echo '  / __ \____  / /_(_)___ ___  (_)___  ____ _/ /_(_)___  ____     ( _ )    '
echo ' / / / / __ \/ __/ / __ `__ \/ /_  / / __ `/ __/ / __ \/ __ \   / __ \/|  '
echo '/ /_/ / /_/ / /_/ / / / / / / / / /_/ /_/ / /_/ / /_/ / / / /  / /_/  <   '
echo '\____/ .___/\__/_/_/ /_/ /_/_/ /___/\__,_/\__/_/\____/_/ /_/   \____/\/   '
echo '    /_/       __            __        __  _                               '
echo '   /   | ____/ /___ _____  / /_____ _/ /_(_)___  ____                     '
echo '  / /| |/ __  / __ `/ __ \/ __/ __ `/ __/ / __ \/ __ \                    '
echo ' / ___ / /_/ / /_/ / /_/ / /_/ /_/ / /_/ / /_/ / / / /                    '
echo '/_/  |_\__,_/\__,_/ .___/\__/\__,_/\__/_/\____/_/ /_/                     '
echo '  ______         /_/____                                                  '
echo ' /_  __/___  ____  / / /_  ____  _  __                                    '
echo '  / / / __ \/ __ \/ / __ \/ __ \| |/_/                                    '
echo ' / / / /_/ / /_/ / / /_/ / /_/ />  <                                      '
echo '/_/  \____/\____/_/_.___/\____/_/|_|                                      '
echo '                                                                          '
echo ' '

# Set goat top directory
setenv LAST_COMMAND `echo $_`
if (`echo ${LAST_COMMAND}` == "") then
  setenv GOATTOP $PWD
else
  setenv SETUP_FILE `echo ${LAST_COMMAND} | cut -d " " -f 2`
  setenv REAL_FILE `eval echo ${SETUP_FILE}`
  setenv SETUP_PATH `dirname ${REAL_FILE}`
  setenv GOATTOP `cd ${SETUP_PATH}; pwd -L`
endif

echo 'Running at goat top directory: '  ${GOATTOP}

# Add scripts to the path variable
setenv GOAT_SCRIPTPATHS ${GOATTOP}/scripts
setenv GOAT_EXECUTABLES ${GOATTOP}/executables
setenv GOAT_VISUALIZATION ${GOATTOP}/scripts
setenv PATH "${GOAT_SCRIPTPATHS}:${GOAT_EXECUTABLES}:${PATH}"

# Set HOST_NAME and COMPILER, which will determine setup files to be used
#------------------------------------------------------------------------
setenv COMPILER gfortran # only gfortran compiler supported (for now)

if (`uname` != "Darwin") then   # Assuming to work on some HPC cluster or on a local Linux device
  if ( $?GOAT_HOST_NAME_FORCE ) then
    setenv HOST_NAME $GOAT_HOST_NAME_FORCE
    echo "Running at $HOST_NAME (set by GOAT_HOST_NAME_FORCE)"
  else if (-s ${GOATTOP}/SETUP/setup.csh.HOST_NAME.local) then
    echo Loading SETUP/setup.csh.HOST_NAME.local.
    source ${GOATTOP}/SETUP/setup.csh.HOST_NAME.local
  else
    if (-s ${GOATTOP}/whereami) then
      set iamat=`${GOATTOP}/whereami|tail -1`
      echo Running at $iamat.
    else
      set iamat="UNKNOWN"
    endif
    switch ($iamat)
    case "*UNKNOWN":
      setenv HOST_NAME LINUX
      breaksw
    case "*LOCAL":
      setenv HOST_NAME LINUX
      breaksw
    default:
      setenv HOST_NAME ${iamat}
    endsw
  endif
else   # Using MacOS, so assuming to work on a local device
  setenv SYSNAME `uname`_`arch`
  setenv ARCH `arch`
  echo Running on MacOS, architecture ${ARCH}.
  setenv HOST_NAME 'DARWIN'
endif

# Load 
set setup=${GOATTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}

# Load default openmp settings
source ${GOATTOP}/SETUP/openmp

# Setup files for combination of HOST_NAME and COMPILER, + local modifications if present
if (-s ${GOATTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}) then
  echo Loading SETUP/setup.csh.${HOST_NAME}.${COMPILER}.
  source ${GOATTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}
else
  echo File SETUP/setup.csh.${HOST_NAME}.${COMPILER} not found!
endif
if (-s ${GOATTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local) then
  echo Loading SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local.
  source ${GOATTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local
endif

# Set some aliases 
alias gtop "cd ${GOATTOP}"
alias pgdinput "python3 ${GOAT_VISUALIZATION}/VisualizeGDInput.py"
alias pgdoutput "python3 ${GOAT_VISUALIZATION}/VisualizeGDOutput.py"
alias mg "python3 ${GOAT_VISUALIZATION}/MonitorGrid.py"
alias psoinput "python3 ${GOAT_VISUALIZATION}/VisualizeShapeOptInput.py"
alias psooutput "python3 ${GOAT_VISUALIZATION}/VisualizeShapeOptOutput.py"
alias mgv "python3 ${GOAT_VISUALIZATION}/MonitorGridAndVessel.py"
alias pgginput "python3 ${GOAT_VISUALIZATION}/VisualizeGGInput.py"
alias pggoutput "python3 ${GOAT_VISUALIZATION}/VisualizeGGOutput.py"
alias pspp "python3 ${GOAT_VISUALIZATION}/SolpsPostProcess.py"
alias pgainput "python3 ${GOAT_VISUALIZATION}/VisualizeGAinput.py"
alias pgaoutput "python3 ${GOAT_VISUALIZATION}/VisualizeGAoutput.py"
alias pgaerror "python3 ${GOAT_VISUALIZATION}/VisualizeGAerror.py"

# List loaded modules
#module list
