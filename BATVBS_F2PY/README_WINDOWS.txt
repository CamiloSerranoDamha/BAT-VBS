******************************************************************************************************************
CONTENT:

Fortran_source_code folder: 
*.f90 files                  : BAT + VBS model Fortran subroutines
f2py_commands_batvbs.bat file: Batch file that contains the command lines to produce a Python 
			       extension module from Fortran source files
Python_examples folder:
*.py files                   : Python scripts that show examples of how to use the BAT + VBS model functions
EXAMPLES_INFO_PYTHON.docx    : Word document that describes the example python files
******************************************************************************************************************
REQUIREMENTS:

* IMPORTANT: Do not mix 32- and 64-bit software. *
If you are using a 64-bit (32-bit) machine, install 64-bit (32-bit) version for all 
of the following:

- Fortran and C compilers (we used gcc 8.3.0): https://strawberryperl.com
- Python (we used Python version 3.9.7, Anaconda Prompt (miniconda3)): 
  https://docs.conda.io/en/latest/miniconda.html
- NumPy package (we used NumPy version 1.21.5, installed with conda): https://anaconda.org/anaconda/numpy 
*******************************************************************************************************************
LIST OF STEPS TO BUILD THE PYTHON EXTENSION MODULE 

1) Produce .o and .mod files, and Python extension module (.pyd)
using F2PY:
Open an Anaconda Prompt (anaconda3 or miniconda3)/activate a conda environment, navigate to the Fortran_source_code directory 
and type the following command: 

f2py_commands_batvbs_windows.bat

[If the procedure is succesful, the Python extension module file (.pyd) should be in your working directory.
If you used Python version 3.x , the name of the extension module file created should contain the reference cp3x.
Make sure to use Python 3.x to import the Python extension module succesfully.]

To run the provided example program from Python, copy the "batvbs... .pyd" file to the folder containing 
the Python examples. Then, execute the Python example like so:  python .\example_bat.py
*******************************************************************************************************************




