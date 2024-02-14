# Binary Activity Thermodynamics (BAT-VBS) model
Fortran version of the Binary Activity Thermodynammics model (BAT-VBS) model (https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model).

If you use the BAT-VBS model, please cite: 

1) Serrano Damha, C., Cummings, B. E., Schervish, M., Shiraiwa, M., Waring, M. S., & Zuend, A. (2024). Capturing the relative-humidity-sensitive gas–particle partitioning of organic aerosols in a 2D volatility basis set. Geophysical Research Letters, 51, e2023GL106095. https://doi.org/10.1029/2023GL106095

2) This repository

3) Gorkowski, K., Preston, T. C., and Zuend, A.: Relative-humidity-dependent organic aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem. Phys., 19, 13383–13407, https://doi.org/10.5194/acp-19-13383-2019, 2019.


CONTENT:

Fortran_source_code folder: 
1) *.f90 files:
BAT and BAT-VBS models Fortran subroutines (source code)

2) EXAMPLE_PROGRAM.f90:
Standalone Fortran program that illustrates how to run the BAT-VBS model

3) BATVBS_F2PY:
Folder containing the procedure to create a BAT-VBS model extension module that can be imported in Python

