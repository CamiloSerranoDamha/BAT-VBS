#****************************************************************************************
#*   example_bat.py                                                                     *
#*                                                                                      *
#*   :: Purpose ::                                                                      *
#*   Script that shows how to use the bat function.                                     *
#*                                                                                      *
#*   bat function                                                                       *
#*   Input  : a binary mixture composed of 1 organic species and its associated         *
#*           water.                                                                     * 
#*                                                                                      *
#*   Output : activity cofficient of water and of the organic species assuming both are *
#*            present in the same single liquid phase (unlike more sophisticated BAT    *
#*            calculations with phase separation consideration).                        *
#*                                                                                      *
#*   :: Author & Copyright ::                                                           *
#*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
#*   McGill University, Montreal, Quebec (2021),                                        *
#*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
#*                                                                                      *
#*   -> created: 15-04-2021                                                             *
#*   -> latest changes: 04-05-2022                                                      *
#!***************************************************************************************

# Import NumPy to specify data types
import numpy as np
# Import the functions of the BAT model
import batvbs

# define the organic species (1 species) that is in the binary system (org+water)
M_org = np.float64(200.05) # Molar mass (g/mol) of the organic species (1 species) in the binary mixture (org+water)
O2C_org = np.float64(0.44) # O/C elemental ratio of the organic species (1 species) in the binary mixture (org+water)
H2C_org = np.float64(1.61) # H/C elemental ratio of the organic species (1 species) in the binary mixture (org+water)
                           # set H/C ratio to a negative real value when it is not known
N2C_org = np.float64(-1.0) # N/C elemental ratio of the organic species (1 species) in the binary mixture (org+water)
                           # set N/C ratio to a negative real value when it is not known
x_org = np.float64(0.75)   # Mole fraction of the organic species (1 species) in the binary mixture (org+water)

Group_org = 3  # Type of oxygen-bearing functional group that represents the organic species (1 species) in the binary mixture (org+water)
# BAT Functional Group Options		
# 1: hydroxyl		
# 2: carboxyl		
# 3: ketone		
# 4: hydroperoxide		
# 5: ether		
# 6: ester		
# 7: hydroperoxideSOA		
# 8: SOA chemicals		
# 9: PEG		

# water properties
M_water = np.float64(18.01528) # Molar mass (g/mol) of water

# calculate the activity cofficient of water and organic species: y_water, y_org
y_water, y_org  = batvbs.bat(M_org, O2C_org, H2C_org, N2C_org, Group_org, x_org) 

# calculate the activity of water and of the organic species (1 species) in the binary mixture (org+water): a_water, a_org
a_water = y_water * (np.float64(1.0) - x_org)   # a_i = y_i * x_i
a_org   = y_org * (x_org)
    
# calculate the mass fraction of water and of the organic species (1 species) in the binary mixture (org+water): w_water, w_org
w_water = ((np.float64(1.0) - x_org) * M_water)/((np.float64(1.0) - x_org) * M_water + x_org * M_org)     # w_i = (x_i * M_i)/(mean_M)
w_org   =  (x_org * M_org)/((np.float64(1.0) - x_org) * M_water + x_org * M_org)

# example output
print("Selected predicted properties for given input")
print("Activity of water and of organic compound : ", format(a_water, '13.6E'), format(a_org, '13.6E'))
print("Mass fraction of water and of organic compound : ", format(w_water, '13.6E'), format(w_org, '13.6E'))