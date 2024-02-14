************************************************************************************************************************
CONTENT:

Fortran_source_code folder: 
*.f90 files                  : BAT + VBS model Fortran subroutines
f2py_commands_batvbs.sh file : Batch file that contains the command lines to produce a Python 
			       extension module from Fortran source files
Python_examples folder:
*.py files                   : Python scripts that show examples of how to use the BAT + VBS model functions
EXAMPLES_INFO_PYTHON.docx    : Word document that describes the example python files
************************************************************************************************************************
REQUIREMENTS:

* IMPORTANT: Do not mix 32- and 64-bit software. *

If you are using a 64-bit (32-bit) machine, install 64-bit (32-bit) version for all 

of the following:


- Fortran and C compilers (we used GCC version 9.3.1 (Linux) and Homebrew GCC 10.2.0 (MacOS))

- Python in a conda environment (we used version 3.9.4 (Linux) and version Python 3.9.2 (MacOS)):  

- NumPy package (we used version 1.20.1 (Linux and MacOS)) installed with conda.

*************************************************************************************************************************
LIST OF STEPS TO BUILD THE PYTHON EXTENSION MODULE USING A CONDA ENVIRONMENT 

1) Make .sh file executable:
Activate a conda environment, navigate to the Fortran_source_code directory and type the following command: 

chmod +x f2py_commands_batvbs_linux_mac.sh


2) Produce .o and .mod files, and Python extension module (.pyd) using F2PY:
Activate a conda environment, navigate to the Fortran_source_code directory and type the following command: 

./f2py_commands_batvbs_linux_mac.sh


[If the procedure is succesful, the Python extension module file (.pyd) should be in your working directory.
If you used Python version 3.x , the name of the extension module file created should contain the reference cp3x.
Make sure to use Python 3.x to import the Python extension module succesfully.]

To run the provided example program from Python, copy the "batvbs... .pyd" file to the folder containing 
the Python examples. Then, execute the Python example like so:  python .\example_bat.py
*************************************************************************************************************************




