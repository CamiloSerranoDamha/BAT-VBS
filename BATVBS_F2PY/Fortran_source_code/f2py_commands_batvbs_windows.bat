@ECHO OFF
SETLOCAL

:: * Creates .o and .mod files
gfortran  -O3 -fPIC -ffree-line-length-none -c  PRECISION_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  DATATYPE_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  MINPACK_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  NN_DATA_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  OBJECTIVE_FUNCTION_VARIABLES_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  TOOLS_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  OBJECTIVE_FUNCTION_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  BAT_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  NN_MOD.f90
gfortran  -O3 -fPIC -ffree-line-length-none -c  BATVBS_MOD.f90

::  ****************************************************************
:: *                                                                *
:: * This is a Windows batch file that wraps the Fortran subroutines* 
:: * of the Binary Activity Thermodynamics (BAT) model              *
:: * (https://doi.org/10.5194/acp-19-13383-2019) to create a Python *
:: * extension module using Fortran to Python Interface Generator   *
:: * (F2PY).                                                        *
:: *                                    (script by Camilo Serrano)  *
::  ****************************************************************  

:: ===== Assign variables for ease of use with commands below =======

:: ================================================================== 
:: The following line specifies the type of the Fortran compiler and the 
:: type of the C compiler that F2PY uses.
:: The following combination of compilers works when installing 
:: Strawbery Perl(https://strawberryperl.com/).

SET _compiler=--fcompiler=gnu95 --compiler=mingw32
:: ================================================================== 


:: ================================================================== 
:: The following line specifies Fortran source code that contains
:: the Python BAT+VBS functions that will be included in the extension
:: module.

SET _f90source=BATVBS_F2PY.f90
:: ================================================================== 


:: ================================================================== 
:: The following line specifies the name of the Python module that will be
:: created.

SET _pythonmodule=batvbs
:: ================================================================== 

:: The following line generates Python extension module (.pyd file):
python -m numpy.f2py -c %_compiler% --f90flags='-Wno-tabs' --f90flags='-O3' *.o %_f90source% -m %_pythonmodule% 

PAUSE

