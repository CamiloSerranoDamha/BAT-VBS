#****************************************************************************************
#*   example_bat_inversion.py                                                           *
#*                                                                                      *
#*   :: Purpose ::                                                                      *
#*   Script that shows how to use the bat_inversion function                            *
#*                                                                                      *
#*   bat_inversion function                                                             *
#*   Input : a mixture of multiple organic species and water.                           *
#*   Output:                                                                            *
#*   1) activity coefficient of each organic species j in liquid phases                 *
#*      alpha (water-rich) and beta (organic-rich),                                     *
#*   2) mass fraction of each organic species j in liquid phases                        *
#*      alpha (water-rich) and beta (organic-rich),                                     *
#*   3) total (organics + water) mass concentration (ug/m3) in liquid phases            * 
#*      alpha (water-rich) and beta (organic-rich).                                     *
#*                                                                                      *
#*   :: Author & Copyright ::                                                           *
#*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
#*   McGill University, Montreal, Quebec (2021),                                        *
#*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
#*                                                                                      *
#*   -> created: 15-04-2021                                                             *
#*   -> latest changes: 04-05-2022                                                      *
#****************************************************************************************

# Import NumpY to specify date types
import numpy as np
# Import the functions of the BAT model
import batvbs

# define the system/simulation properties (a system of 10 organic species and water in this example)
a_water = np.float64(0.4000) # Target water activity (RH) that the mole fraction of water (and organics) in the mixture must match
M_water = np.float64(18.01528) # Molar mass (g/mol) of water

# The first element corresponds to the first organic species in the system, the second element corresponds to the second organic species j in the system, etc
M_org_j = np.array([200.0, 188.0, 216.0, 368.0, 368.0, 204.0, 195.0, 368.0, 158.0, 206.0], np.float64) # Molar mass of each organic species j in the system (g/mol)
O2C_org_j = np.array([0.4, 0.444, 0.5, 0.368, 0.368, 0.556, 0.857, 0.368, 0.375, 0.75], np.float64) # Elemental O/C ratio of each organic species j in the system
H2C_org_j = np.array([1.6 , 1.78, 1.6, 1.47, 1.56, 1.78, 1.75, 1.56, 1.75, 1.75], np.float64) # Elemental H/C ratio of each organic species j in the system
                                                                                              # H/C is optional; set H/C to a negative real number if it is not known
N2C_org_j = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0], np.float64)# Elemental N/C ratio of each organic species j in the system
                                                                                              # N/C is optional; set N/C to a negative real number if it is not known
C_OA_org_j = np.array([4.79E+00 , 1.98E+00 , 0.86E+00 , 2.07E+00 , 0.72E+00 , 0.519E+00, 0.366E+00, 0.82E+00 , 0.199E+00, 0.213E+00], np.float64) # Total "dry" mass concentration of each organic species j in the particle (ug/m3)
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

# calculate the activity cofficient of each organic species j in liquid phases alpha and beta: y_org_j_alpha, y_org_j_beta, 
# calculate the mass fraction of each organic species j in liquid phases alpha and beta: w_org_j_alpha, w_org_j_beta,
# calculate the total (organic + water) mass concentration in liquid phases alpha and beta (ug/m3): C_OA_water_org_alpha, C_OA_water_org_beta                                 
y_org_j_alpha, y_org_j_beta, w_org_j_alpha, w_org_j_beta, C_OA_water_org_alpha, C_OA_water_org_beta = batvbs.bat_inversion(M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, Group_org_j, C_OA_org_j, a_water) 

# calculate the mass concentration of each organic species j in liquid phases alpha and beta (ug/m3): C_OA_org_j_alpha, C_OA_org_j_beta
C_OA_org_j_alpha = C_OA_water_org_alpha * w_org_j_alpha
C_OA_org_j_beta = C_OA_water_org_beta * w_org_j_beta

# calculate the total organic mass concentration in liquid phases alpha and beta (ug/m3): C_OA_org_alpha, C_OA_org_beta
C_OA_org_alpha = np.sum(C_OA_org_j_alpha)
C_OA_org_beta = np.sum(C_OA_org_j_beta) 

# calculate the mass fraction of water in liquid phases alpha and beta: w_water_alpha, w_water_beta
w_water_alpha = np.float64(1.0) - np.sum(w_org_j_alpha)
w_water_beta = np.float64(1.0) - np.sum(w_org_j_beta)

# calculate the water mass concentration in liquid phases alpha and beta (ug/m3): C_OA_water_alpha, C_OA_water_beta
C_OA_water_alpha = C_OA_water_org_alpha * w_water_alpha
C_OA_water_beta = C_OA_water_org_beta * w_water_beta

# calculate the total (organic + water; liquid phase alpha + beta) mass concentration in the aqueous organic particle (ug/m3): C_OA_water_org_alpha_beta
C_OA_water_org_alpha_beta = C_OA_water_org_alpha + C_OA_water_org_beta
# also, C_OA_water_org_alpha_beta = (C_OA_org_alpha + C_OA_water_alpha) + (C_OA_org_beta + C_OA_water_beta)

# calculate the fractional liquid–liquid partitioning of each organic species j to liquid phases alpha and beta: q_org_j_alpha, q_org_j_beta
# in the case of a liquid–liquid equilibrium, the relative phase preferences are described by q_org_j_alpha and q_org_j_beta
q_org_j_alpha = C_OA_org_j_alpha/C_OA_org_j
q_org_j_beta = np.float64(1.0) - q_org_j_alpha

# calculate the instantaneous effective saturation concentration C* according to equation 6 in https://doi.org/10.5194/acp-19-13383-2019: C_star_org_j
C_sat_org_j = np.array([5.74E+03, 3.27E+02, 1.67E+02, 2.79E-06, 1.05E+02, 2.13E+00, 7.19E-01, 3.64E-06, 1.16E+03, 3.02E-02], np.float64) # Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
# via liquid phase alpha
C_star_org_j_alpha = C_sat_org_j * y_org_j_alpha * q_org_j_alpha * C_OA_water_org_alpha_beta/(M_org_j * (np.sum(C_OA_org_j_alpha/(M_org_j)) + (C_OA_water_alpha/M_water)))
# via liquid phase beta
C_star_org_j_beta = C_sat_org_j * y_org_j_beta * q_org_j_beta * C_OA_water_org_alpha_beta/(M_org_j * (np.sum(C_OA_org_j_beta/(M_org_j)) + (C_OA_water_beta/M_water)))   
# weighted mean using q_org_j_alpha values as weights
C_star_org_j = (C_star_org_j_alpha * q_org_j_alpha) + (C_star_org_j_beta * q_org_j_beta)

# calculate the hygroscopicity parameter kappa of the organic aerosol only from BAT water uptake calculations (before VBS): kappaHGF
# estimate the density of each organic species j (g/cm3): density_org_j
# set H2C_org_j/N2C_org_j to a negative float number if its value is not known (these are the only optional inputs to the function below)
density_org_j = batvbs.density_org_estimate(M_org_j, O2C_org_j, H2C_org_j, N2C_org_j)
# density of water (g/cm3): density_water
density_water = np.float64(0.997)
# calculate the cumulative contribution of organic component volumes at a given RH (ug cm3 m-3 g-1): V_org
V_org = np.sum((C_OA_org_j_alpha + C_OA_org_j_beta)/density_org_j)
# calculate the cumulative contribution of water volume at given RH (ug cm3 m-3 g-1): V_water
V_water = (C_OA_water_alpha + C_OA_water_beta)/density_water
# calculate the hygroscopicity parameter kappa of OA:
kappaHGF = ((np.float64(1.0)/a_water)-np.float64(1.0))*(V_water/V_org)

# example output
print("Selected predicted properties for given input at RH = ", format(a_water, '13.6E'))
print("Mass concentration of org and water in liquid phase alpha (ug/m3): ", format(C_OA_water_org_alpha, '13.6E'))
print("Mass concentration of org and water in liquid phase beta (ug/m3): ", format(C_OA_water_org_beta, '13.6E'))
print("Mass concentration of org and water in the particle (phase alpha + beta) (ug/m3): ", format(C_OA_water_org_alpha_beta, '13.6E')) 
np.set_printoptions(formatter={'float_kind':'{:13.6E}'.format})   
print("Instantaneous effective saturation concentration, C*_j, of each organic species j (ug/m3): ", C_star_org_j)
print("Hygroscopicity parameter kappa of the organic aerosol: ", format(kappaHGF, '13.6E'))
