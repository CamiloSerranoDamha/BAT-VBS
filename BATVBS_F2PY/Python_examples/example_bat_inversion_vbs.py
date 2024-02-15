#*****************************************************************************************
#*   example_bat_inversion_vbs.py                                                        *
#*                                                                                       *
#*   :: Purpose ::                                                                       *
#*   Script used to test the Python functions created through F2PY.                      *
#*                                                                                       *
#*   bat_inversion_vbs function                                                          *
#!*   Input system: A mixture of multiple organic species and water.                     *
#!*   Output:                                                                            *
#!*   1) equilibrium mass concentration (ug/m3) of each organic species j in liquid      *
#!*      phases alpha (water-rich) and beta (organic-rich) that results from BAT+VBS     *
#!*      calculations (ug/m3),                                                           *
#!*   2) equilibrium water mass concentration (ug/m3) in liquid phases alpha (water-rich)* 
#!*      and beta (organic-rich) that results from BAT+VBS calculations.                 *
#!*                                                                                      *                                                                    *
#*   :: Author & Copyright ::                                                            *
#*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                       *
#*   McGill University, Montreal, Quebec (2021),                                         *
#*   Dept. Atmospheric and Oceanic Sciences, McGill University                           *
#*                                                                                       *
#*   -> created: 15-04-2021                                                              *
#*   -> latest changes: 04-05-2022                                                       *
#*****************************************************************************************

# Import NumpY to specify date types
import numpy as np
# Import the functions of the BAT model
import batvbs

# define the system/simulation properties (a system of 10 organic species and water in this example)
a_water = np.float64(0.4000) # Target water activity (RH) that the mole fraction of water (and organics) in the mixture must match

# The first element corresponds to the first organic species in the system, the second element corresponds to the second organic species j in the system, etc
M_org_j = np.array([200.0, 188.0, 216.0, 368.0, 368.0, 204.0, 195.0, 368.0, 158.0, 206.0], np.float64) # Molar mass of each organic species j in the system (g/mol)
O2C_org_j = np.array([0.4, 0.444, 0.5, 0.368, 0.368, 0.556, 0.857, 0.368, 0.375, 0.75], np.float64) # Elemental O/C ratio of each organic species j in the system
H2C_org_j = np.array([1.6 , 1.78, 1.6, 1.47, 1.56, 1.78, 1.75, 1.56, 1.75, 1.75], np.float64) # Elemental H/C ratio of each organic species j in the system
                                                                                              # H/C is optional; set H/C to a negative real number if it is not known
N2C_org_j = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0], np.float64)# Elemental N/C ratio of each organic species j in the system
                                                                                              # N/C is optional; set N/C to a negative real number if it is not known
C_gas_particle_org_j = np.array([8.79, 3.98, 1.13, 4.07, 1.02, 0.919, 0.766, 1.02, 0.399, 0.313], np.float64) # Mass concentration of each organic species j in gas phase + particle phase (ug/m3)
C_sat_org_j = np.array([5.74E+03, 3.27E+02, 1.67E+02, 2.79E-06, 1.05E+02, 2.13E+00, 7.19E-01, 3.64E-06, 1.16E+03, 3.02E-02], np.float64) # Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
Group_org_j = [3,3,3,3,3,3,3,3,3,3] # Type of oxygen-bearing functional group of each organic species j in the system
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

# At a given RH:
# calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_org_j_alpha, C_OA_org_j_alpha,
# calculate the equilibrium water mass concentration in liquid phases alpha and beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_water_alpha, C_OA_water_beta
C_OA_org_j_alpha, C_OA_org_j_beta, C_OA_water_alpha, C_OA_water_beta = batvbs.bat_inversion_vbs(M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_org_j, C_gas_particle_org_j, Group_org_j, a_water)

# calculate the total (water + organic) equilibrium mass concentration in liquid phases alpha and beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_water_org_alpha, C_OA_water_org_beta
C_OA_water_org_alpha = np.sum(C_OA_org_j_alpha) + C_OA_water_alpha
C_OA_water_org_beta = np.sum(C_OA_org_j_beta) + C_OA_water_beta

# calcualte the total (organics + water; liquid phase alpha + beta) mass concentration in the aqueous organic particle that results from VBS calculations (ug/m3): C_OA_water_org_alpha_beta
C_OA_water_org_alpha_beta = C_OA_water_org_alpha + C_OA_water_org_beta

# calculate the hygroscopicity parameter kappa of OA: kappaHGF
# this is an estimation of kappa of the OA based on the outputs of the coupled BAT+VBS model
# estimate the density of each organic species j (g/cm3): density_org_j
# set H2C_org_j/N2C_org_j to a negative float number if its value is not known (these are the only optional inputs to the function below)
density_org_j = batvbs.density_org_estimate(M_org_j, O2C_org_j, H2C_org_j, N2C_org_j)
# density of water (g/cm3): density_water
density_water = np.float64(0.997)
# calculate the cumulative contribution of organic component volumes at a given RH (ug cm3 m-3 g-1): V_org
V_org = np.sum((C_OA_org_j_alpha + C_OA_org_j_beta)/density_org_j)

# calculate the cumulative contribution of water volume at a given RH (ug cm3 m-3 g-1): V_water
V_water = (C_OA_water_alpha + C_OA_water_beta)/density_water
# calculate the hygroscopicity parameter kappa of OA:
kappaHGF = ((np.float64(1.0)/a_water)-np.float64(1.0))*(V_water/V_org)

# example output
print("Selected predicted properties for given input at RH = ", format(a_water, '13.6E'))
print("Equilibrium mass concentration of org + water in liquid phase alpha (ug/m3): ", format(C_OA_water_alpha + np.sum(C_OA_org_j_alpha), '13.6E'))
print("Equilibrium mass concentration of org + water in liquid phase beta (ug/m3): ", format(C_OA_water_beta + np.sum(C_OA_org_j_beta), '13.6E'))
print("Equilibrium mass concentration of org + water in the particle (ug/m3): ", format(C_OA_water_org_alpha_beta, '13.6E'))
print("Equilibrium mass concentration of water in the particle (ug/m3): ", format(C_OA_water_alpha + C_OA_water_beta, '13.6E'))
print("Equilibrium mass concentration of organics in the particle (ug/m3): ", format(np.sum(C_OA_org_j_alpha) + np.sum(C_OA_org_j_beta), '13.6E'))
print("Hygroscopicity parameter kappa of the organic aerosol: ", format(kappaHGF, '13.6E'))
