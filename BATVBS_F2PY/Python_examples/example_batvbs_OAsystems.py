# -*- coding: utf-8 -*-

#****************************************************************************************#
#*   BATVBS_testcase.py                                                                 *#
#*                                                                                      *#
#*   :: Purpose ::                                                                      *#
#*   Script used to estimate the equilibrium organic and org. water mass                *#
#*   concentrations in the aqueous org. particle using the coupled BAT+VBS model.       *#
#*                                                                                      *#
#*   bat_inversion_vbs function                                                         *#
#*   Input system: A mixture of multiple organic species j and water.                   *#
#*   Output:                                                                            *#
#*   1) equilibrium mass concentration (ug/m3) of each organic species j in liquid      *#
#*      phases alpha (water-rich) and beta (organic-rich) that results from BAT+VBS     *#
#*      calculations (ug/m3),                                                           *#
#*   2) equilibrium water mass concentration (ug/m3) in liquid phases alpha (water-rich)*#
#*      and beta (organic-rich) that results from BAT+VBS calculations.                 *#
#*                                                                                      *#
#*   :: Author & Copyright ::                                                           *#
#*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *#
#*   McGill University, Montreal, Quebec (2021),                                        *#
#*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *#
#*                                                                                      *#
#*   -> created: 15-04-2021                                                             *#
#*   -> latest changes: 02-02-2022                                                      *#
#****************************************************************************************#

# Import NumpY to specify date types
import numpy as np
# Import the functions of the BAT model
import batvbs
# Import pandas to read Excel files
import pandas as pd
# Import matplotlib.pyplot plotting library
import matplotlib.pyplot as plt

# Path of the indoor OA data file
loc = 'Input_OA/VBS-BAT_prelim-inputs.xlsx'

# Read all the sheets of the Workbook
df = pd.read_excel(loc, index_col=0, sheet_name=None, na_filter=False)

#################################
### STEP 1: DEFINE THE SYSTEM ###
#################################

## Shared system properties (Same for HOA, OOA and HOA+OOA systems) ##

# The coupled BAT+VBS model will be run at different RH values (from 0.0 to 0.95 with a step of 0.01)
# Target water activity (RH) that the mole fraction of water (and organics) in the mixture must match
a_water = np.arange(0.00, 0.96, 0.01, np.float64)

# Define possible values of O/C ratio and pure component saturation mass concentration at the temperature of interest (ug/m3)
O2C_possible_values = np.array(
    [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], np.float64)
C_star_Tref_possible_values = np.array(
    [10.0**-7, 10.0**-6, 10.0**-5, 10.0**-4, 10.0**-3, 10.0**-2, 10.0**-1, 10.0**0, 10.0**1, 10.0**2, 10.0**3], np.float64)

# Molar mass data (g/mol) of each organic species j in the system
M_org_j = ((df['molar mass']).to_numpy(np.float64)).flatten()

# Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
# Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
C_star_Tref_org_j = np.array(np.resize(C_star_Tref_possible_values, np.size(
    O2C_possible_values)*np.size(O2C_possible_values)), np.float64)

#T_ref = # K
#enthalpy_vap = # kJ/mol
#R = # kJ/(mol K)
# Apply Clausius-Clapeyron equation if needed here
#C_star_T_org_j = C_star_Tref_org_j*(T_ref/T)*exp(-(enthalpy_vap/R)*(1/T_ref-1/T))
# Otherwise:
C_star_T_org_j = C_star_Tref_org_j

# Elemental O/C ratio of each organic species j in the system
O2C_org_j = np.array(np.repeat(O2C_possible_values,
                     np.size(C_star_Tref_possible_values)))

# Elemental H/C ratio of each organic species j in the system
# set H/C to a negative real number (not known)
H2C_org_j = np.ones(np.size(M_org_j), np.float64) * -1.0

# Elemental N/C ratio of each organic species j in the system
# set N/C to a negative real number (not known)
N2C_org_j = np.ones(np.size(M_org_j), np.float64) * -1.0

# Type of oxygen-bearing functional group of each organic species j in the system
Group_org_j = np.repeat([3], np.size(M_org_j))
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

# Estimate the density of each organic species j (g/cm3): density_org_j
# set H2C_org_j/N2C_org_j to a negative float number if its value is not known (these are the only optional inputs to the function below)
density_org_j = batvbs.density_org_estimate(M_org_j, O2C_org_j, H2C_org_j, N2C_org_j)
# density of water (g/cm3): density_water
density_water = np.float64(0.997)

## Specific system properties (different for HOA, OOA and HOA+OOA systems) ##
# Total (gas + particle) mass concentration data (ug/m3) of each organic species j in the HOA system
C_gas_particle_org_j_HOA = ((df['HOA']).to_numpy(np.float64)).flatten()
#data_HOA = np.nonzero(C_gas_particle_org_j_HOA)
# Total (gas + particle) mass concentration data (ug/m3) of each organic species j in the OOA system
C_gas_particle_org_j_OOA = ((df['OOA']).to_numpy(np.float64)).flatten()
#data_OOA = np.nonzero(C_gas_particle_org_j_OOA)
# Total (gas + particle) mass concentration data (ug/m3) of each organic species j in the HOA+OOA system
C_gas_particle_org_j_HOA_OOA = ((df['HOA+OOA']).to_numpy(np.float64)).flatten()
#data_HOA_OOA = np.nonzero(C_gas_particle_org_j_HOA_OOA)

########################################
### STEP 2: SIMULATION USING BAT+VBS ###
########################################

# Create output arrays
# Equilibrium organic mass concentration of each org. species in the particle phase (phases alpha + beta) that results from coupled BAT+VBS calculations (ug/m3)
C_OA_org_j_HOA = np.zeros([np.size(a_water), np.size(M_org_j)], np.float64)
C_OA_org_j_OOA = np.zeros([np.size(a_water), np.size(M_org_j)],np.float64)
C_OA_org_j_HOA_OOA = np.zeros([np.size(a_water), np.size(M_org_j)], np.float64)
# Equilibrium organic mass concentration (sum of all org. species) in the particle phase (phases alpha + beta) that results from coupled BAT+VBS calculations (ug/m3)
C_OA_org_HOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_OOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_HOA_OOA = np.zeros(np.size(a_water), np.float64)
# Equilibrium organic mass concentration (sum of all org. species) in the particle phase alpha and particle phase beta that results from coupled BAT+VBS calculations (ug/m3)
C_OA_org_alpha_HOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_beta_HOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_alpha_OOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_beta_OOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_alpha_HOA_OOA = np.zeros(np.size(a_water), np.float64)
C_OA_org_beta_HOA_OOA = np.zeros(np.size(a_water), np.float64)
# Equilibrium water mass concentration (due to organics) in the particle phase (phases alpha + beta) that results from coupled BAT+VBS calculations (ug/m3)
C_OA_water_HOA = np.zeros(np.size(a_water), np.float64)
C_OA_water_OOA = np.zeros(np.size(a_water), np.float64)
C_OA_water_HOA_OOA = np.zeros(np.size(a_water), np.float64)
# Hygroscopicity parameter kappa of OA, related to the hygroscopic growth factor
kappa_HOA = np.zeros(np.size(a_water), np.float64)
kappa_OOA = np.zeros(np.size(a_water), np.float64)
kappa_OOA = np.zeros(np.size(a_water), np.float64)
kappa_HOA_OOA = np.zeros(np.size(a_water), np.float64)

# Start simulation

# estimate the pure component saturation concentration at the temperature of interest (C_sat_possible_values) 
C_sat_T_org_j_HOA = batvbs.csat_approx(M_org_j, O2C_org_j, C_star_T_org_j, C_gas_particle_org_j_HOA, Group_org_j)
C_sat_T_org_j_OOA = batvbs.csat_approx(M_org_j, O2C_org_j, C_star_T_org_j, C_gas_particle_org_j_OOA, Group_org_j)
C_sat_T_org_j_HOA_OOA = batvbs.csat_approx(M_org_j, O2C_org_j, C_star_T_org_j, C_gas_particle_org_j_HOA_OOA, Group_org_j)
    
for i in range(0, np.size(a_water)):
    
    # calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_org_j_alpha, C_OA_org_j_alpha,
    # calculate the equilibrium water mass concentration in liquid phases alpha and beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_water_alpha, C_OA_water_beta
    C_OA_org_j_alpha_HOA, C_OA_org_j_beta_HOA, C_OA_water_alpha_HOA, C_OA_water_beta_HOA = batvbs.bat_inversion_vbs(
        M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_T_org_j_HOA, 
        C_gas_particle_org_j_HOA, Group_org_j, a_water[i])
    
    C_OA_org_j_alpha_OOA, C_OA_org_j_beta_OOA, C_OA_water_alpha_OOA, C_OA_water_beta_OOA = batvbs.bat_inversion_vbs(
        M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_T_org_j_OOA, 
        C_gas_particle_org_j_OOA, Group_org_j, a_water[i])
    
    C_OA_org_j_alpha_HOA_OOA, C_OA_org_j_beta_HOA_OOA, C_OA_water_alpha_HOA_OOA, C_OA_water_beta_HOA_OOA = batvbs.bat_inversion_vbs(
        M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_T_org_j_HOA_OOA, 
        C_gas_particle_org_j_HOA_OOA, Group_org_j, a_water[i])

    # calculate the equilibrium organic mass concentration (sum of all org. species) in the particle phase (phase alpha, phase beta, phase alpha+beta) that results from coupled BAT+VBS calculations (ug/m3): C_OA_org
    # calculate the equilibrium water mass concentration in liquid phase alpha+beta that results from coupled BAT+VBS calculations (ug/m3): C_OA_water
    C_OA_org_HOA[i] = np.sum(C_OA_org_j_alpha_HOA + C_OA_org_j_beta_HOA)
    C_OA_org_j_HOA[i] = (C_OA_org_j_alpha_HOA + C_OA_org_j_beta_HOA)
    C_OA_org_alpha_HOA[i] = np.sum(C_OA_org_j_alpha_HOA)
    C_OA_org_beta_HOA[i] = np.sum(C_OA_org_j_beta_HOA)
    C_OA_water_HOA[i] = C_OA_water_alpha_HOA + C_OA_water_beta_HOA
    C_OA_org_OOA[i] = np.sum(C_OA_org_j_alpha_OOA + C_OA_org_j_beta_OOA)
    C_OA_org_j_OOA[i] = (C_OA_org_j_alpha_OOA + C_OA_org_j_beta_OOA)
    C_OA_org_alpha_OOA[i] = np.sum(C_OA_org_j_alpha_OOA)
    C_OA_org_beta_OOA[i] = np.sum(C_OA_org_j_beta_OOA)
    C_OA_water_OOA[i] = C_OA_water_alpha_OOA + C_OA_water_beta_OOA
    C_OA_org_HOA_OOA[i] = np.sum(C_OA_org_j_alpha_HOA_OOA + C_OA_org_j_beta_HOA_OOA)
    C_OA_org_j_HOA_OOA[i] = (C_OA_org_j_alpha_HOA_OOA + C_OA_org_j_beta_HOA_OOA)
    C_OA_org_alpha_HOA_OOA[i] = np.sum(C_OA_org_j_alpha_HOA_OOA)
    C_OA_org_beta_HOA_OOA[i] = np.sum(C_OA_org_j_beta_HOA_OOA)
    C_OA_water_HOA_OOA[i] = C_OA_water_alpha_HOA_OOA + C_OA_water_beta_HOA_OOA
    
    # the following volumes will be used to estimate the hygroscopicity parameter kappa of OA (KappaHGF)
    # calculate the cumulative contribution of organic component volumes at different RH values (ug cm3 m-3 g-1): V_org
    V_org_HOA = np.sum((C_OA_org_j_alpha_HOA + C_OA_org_j_beta_HOA)/(density_org_j))
    V_org_OOA = np.sum((C_OA_org_j_alpha_OOA + C_OA_org_j_beta_OOA)/(density_org_j))   
    V_org_HOA_OOA = np.sum((C_OA_org_j_alpha_HOA_OOA + C_OA_org_j_beta_HOA_OOA)/(density_org_j)) 
    # calculate the cumulative contribution of water volume at different RH values (ug cm3 m-3 g-1): V_water
    V_water_HOA = (C_OA_water_alpha_HOA + C_OA_water_beta_HOA)/density_water
    V_water_OOA = (C_OA_water_alpha_OOA + C_OA_water_beta_OOA)/density_water
    V_water_HOA_OOA = (C_OA_water_alpha_HOA_OOA + C_OA_water_beta_HOA_OOA)/density_water
    
    # calculate the hygroscopicity parameter kappa of OA:
    # ** Warming: You will get division by 0 error at RH = 0 **
    kappa_HOA[i] = ((np.float64(1.0)/a_water[i])-np.float64(1.0)) * ((V_water_HOA)/V_org_HOA)
    kappa_OOA[i] = ((np.float64(1.0)/a_water[i])-np.float64(1.0)) * ((V_water_OOA)/V_org_OOA)
    kappa_HOA_OOA[i] = ((np.float64(1.0)/a_water[i])-np.float64(1.0)) * ((V_water_HOA_OOA)/V_org_HOA_OOA)       



# ########################################
############ STEP 3: PLOTS #############
########################################

# C_OA_org
plt.plot(a_water, C_OA_org_HOA, color='grey', label='HOA', linewidth=1.0)
plt.plot(a_water, C_OA_org_HOA_OOA, color='yellow', label='HOA & OOA', linewidth=1.0)
plt.plot(a_water, C_OA_org_OOA, color='green', label='OOA', linewidth=1.0)
plt.ylim(5.0, 7.5)
plt.ylabel('OA organic mass [$\mathrm{\mu g~m^{-3}}$]')
plt.yticks(np.arange(5.0, 7.75, step=0.25))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper left')
plt.savefig('Output_OA/C_OA_org.pdf', dpi=300)
#plt.show()
plt.close()

plt.plot(a_water, C_OA_org_alpha_HOA, color='blue', label='HOA-alpha', linewidth=1.0)
plt.plot(a_water, C_OA_org_beta_HOA, color='green', label='HOA-beta', linewidth=1.0)
plt.ylim(5.0, 7.5)
plt.ylabel('OA organic mass [$\mathrm{\mu g~m^{-3}}$]')
plt.yticks(np.arange(5.0, 7.75, step=0.25))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper left')
plt.savefig('Output_OA/C_OA_org_HOA.pdf', dpi=300)
#plt.show()
plt.close()

plt.plot(a_water, C_OA_org_alpha_OOA, color='blue', label='OOA-alpha', linewidth=1.0)
plt.plot(a_water, C_OA_org_beta_OOA, color='green', label='OOA-beta', linewidth=1.0)
plt.ylim(5.0, 7.5)
plt.ylabel('OA organic mass [$\mathrm{\mu g~m^{-3}}$]')
plt.yticks(np.arange(5.0, 7.75, step=0.25))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper left')
plt.savefig('Output_OA/C_OA_org_OOA.pdf', dpi=300)
#plt.show()
plt.close()

plt.plot(a_water, C_OA_org_alpha_HOA_OOA, color='blue', label='HOA&OOA-alpha', linewidth=1.0)
plt.plot(a_water, C_OA_org_beta_HOA_OOA, color='green', label='HOA&OOA-beta', linewidth=1.0)
plt.ylim(5.0, 7.5)
plt.ylabel('OA organic mass [$\mathrm{\mu g~m^{-3}}$]')
plt.yticks(np.arange(5.0, 7.75, step=0.25))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper left')
plt.savefig('Output_OA/C_OA_org_HOA_OOA.pdf', dpi=300)
#plt.show()
plt.close()

# C_OA_water
plt.plot(a_water, C_OA_water_HOA, color='grey', label='HOA', linewidth=1.0)
plt.plot(a_water, C_OA_water_HOA_OOA, color='yellow', label='HOA & OOA', linewidth=1.0)
plt.plot(a_water, C_OA_water_OOA, color='green', label='OOA', linewidth=1.0)
plt.ylim(0.0, 6.0)
plt.ylabel('OA water mass [$\mathrm{\mu g~m^{-3}}$]')
plt.yticks(np.arange(0.0, 6.5, step=0.5))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper left')
plt.savefig('Output_OA/C_OA_water.pdf', dpi=300)
#plt.show()
plt.close()

# KappaOA
plt.plot(a_water, kappa_HOA, color='grey', label='HOA', linewidth=1.0)
plt.plot(a_water, kappa_HOA_OOA, color='yellow', label='HOA & OOA', linewidth=1.0)
plt.plot(a_water, kappa_OOA, color='green', label='OOA', linewidth=1.0)
plt.ylim(0.0, 0.15)
plt.ylabel('$\u03BA^{\mathrm{OA}}_{org}$')
plt.yticks(np.arange(0.0, 0.16, step=0.01))
plt.xlim(0, 1.0)
plt.xlabel('$\mathrm{RH}$')
plt.xticks(np.arange(0, 1.1, step=0.1)) 
plt.legend(loc='upper right')
plt.savefig('Output_OA/Kappa_OA.pdf', dpi=300)
#plt.show()
plt.close()

########################################
######### STEP 4: Save Output ##########
########################################

# creating dataframe 
df_HOA = pd.DataFrame({
    'RH':a_water, 
    'C_OA_org (ug/m3)':C_OA_org_HOA, 
    'C_OA_org_alpha (ug/m3)':C_OA_org_alpha_HOA, 
    'C_OA_org_beta (ug/m3)':C_OA_org_beta_HOA,     
    'C_OA_water (ug/m3)': C_OA_water_HOA, 
    'kappa_OA_org': kappa_HOA})

df_OOA = pd.DataFrame({
    'RH':a_water, 
    'C_OA_org (ug/m3)':C_OA_org_OOA, 
    'C_OA_org_alpha (ug/m3)':C_OA_org_alpha_OOA, 
    'C_OA_org_beta (ug/m3)':C_OA_org_beta_OOA,     
    'C_OA_water (ug/m3)': C_OA_water_OOA, 
    'kappa_OA_org': kappa_OOA})

df_HOA_OOA = pd.DataFrame({
    'RH':a_water, 
    'C_OA_org (ug/m3)':C_OA_org_HOA_OOA, 
    'C_OA_org_alpha (ug/m3)':C_OA_org_alpha_HOA_OOA, 
    'C_OA_org_beta (ug/m3)':C_OA_org_beta_HOA_OOA,     
    'C_OA_water (ug/m3)': C_OA_water_HOA_OOA, 
    'kappa_OA_org': kappa_HOA_OOA})

# save to .xlsx file
df_HOA.to_excel('Output_OA/HOA_org_water.xlsx', index=False, sheet_name='HOA')
df_OOA.to_excel('Output_OA/OOA_org_water.xlsx', index=False, sheet_name='OOA')
df_HOA_OOA.to_excel('Output_OA/HOA_OOA_org_water.xlsx', index=False, sheet_name='HOA+OOA')

