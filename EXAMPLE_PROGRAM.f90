!****************************************************************************
!  PROGRAM: MAIN_PROG
!
!  PURPOSE:  TEST BAT+VBS MODEL WITH HOA, OOA and HOA&OOA SYSTEMS
!
!  CREATED: 2022-11-24
!  LAST MODIFICATION: 2022-11-24
!****************************************************************************

PROGRAM MAIN_PROG

    IMPLICIT NONE
    
    !On variable KIND: For compatibility with typical choices in other Fortran programs and f2py,  
    !we use INTEGER(ip) and REAL(rp), where ip = KIND(1.0E0) and rp = KIND(1.0D0), as precision KINDs in this main program (serving as interface). 
    !Within the BAT procedures, a set working real precision (wrp) and similar for integer is used, 
    !which are defined in module PRECISION_MOD (and accessible via BATVBS_MOD).
    
    ! Notation: 
    ! org_j : each individual organic species j
    ! org : total organics 
    ! water: total water
    ! water_org: total water + total organics
    ! alpha: in particle (liquid) phase alpha (water-rich)
    ! beta: in particle (liquid) phase beta (organic-rich)
    ! alpha_beta: in particle (liquid) phase alpha + in particle (liquid) phase beta
    ! dry: variables related to BAT+VBS calculations under dry conditions (RH = activity_water = 0.0)
    
    ! Inputs
    INTEGER, PARAMETER                      :: ip = KIND(1.0E0)                                         ! Precision of INTEGERS
    INTEGER, PARAMETER                      :: rp = KIND(1.0D0)                                         ! Precision of REALS
    INTEGER(ip), PARAMETER					:: N_org = 20                                               ! Number of organic species in the system
	INTEGER(ip), PARAMETER                  :: N_a_water = 1	                                        ! Number of water activity/RH values
    INTEGER(ip), DIMENSION(N_org)           :: Group_org_j                                              ! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(rp), DIMENSION(N_org)              :: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, &              ! Molecular weight (g/mol) and elemental ratios (O/C, H/C, N/C) of each organic species j in the system
                                                                                                        ! * H/C and N/C are optional inputs; set their value to a negative real number when it is not known *
                                                & C_gas_particle_org_j, &                               ! Mass concentration of each organic species j in gas phase + particle phase (ug/m3)
                                                & C_star_org_j, &										! Effective saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
	                                            & C_sat_org_j                                           ! Pure-component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
    REAL(rp), DIMENSION(N_a_water)			:: activity_water                                           ! Target relative humidity (gas-phase water activity) that the mole fraction of water (and organics) in the mixture must match (scale is from 0 to 1)

    ! Local variables
    ! mass concentrations at a given water activity (RH)
    REAL(rp), DIMENSION(N_a_water)			:: C_PM_water_alpha, C_PM_water_beta                        ! Equilibrium water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at any RH > 0.0 (ug/m3)
    REAL(rp), DIMENSION(N_a_water, N_org)	:: C_PM_org_j_alpha, C_PM_org_j_beta                        ! Equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at any RH > 0.0 (ug/m3)

    ! set relative humidity values
    activity_water = [0.38_rp]

    ! the first element corresponds to the first organic species in the system, the second element corresponds to the second organic species in the system, etc
    M_org_j = [528.482_rp, 495.767_rp, 463.052_rp, 430.336_rp, 397.621_rp, 364.906_rp, 332.190_rp, 299.475_rp, 463.052_rp, 234.044_rp, 530.408_rp, &
		& 497.573_rp, 464.739_rp, 431.904_rp, 399.069_rp, 366.235_rp, 333.400_rp, 300.566_rp, 464.739_rp, 234.897_rp]
    O2C_org_j = [2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 2.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, &
		& 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp, 1.0E-01_rp]  ! O/C is not optional
    H2C_org_j = [-1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp,-1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, &
		& -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp]	! H/C is optional; set H/C to a negative real number if it is not known
    N2C_org_j = [-1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp,-1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, &
		& -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp, -1.0E+00_rp]        ! N/C is optional; set N/C to a negative real number if it is not known 
    C_star_org_j = [1.0E-06_rp, 1.0E-05_rp, 1.0E-04_rp, 1.0E-03_rp, 1.0E-02_rp, 1.0E-01_rp, 1.0E+00_rp, 1.0E+01_rp, 1.0E+02_rp, 1.0E+03_rp, 1.0E-06_rp, 1.0E-05_rp, 1.0E-04_rp, 1.0E-03_rp, &
		& 1.0E-02_rp, 1.0E-01_rp, 1.0E+00_rp, 1.0E+01_rp, 1.0E+02_rp, 1.0E+03_rp]
    C_gas_particle_org_j  = [0.0886549E+00_rp, 0.0940607E+00_rp, 0.10271E+00_rp, 0.12109E+00_rp, 0.15893E+00_rp, 0.235692E+00_rp, 0.389217E+00_rp, 0.697347E+00_rp, 1.31901E+00_rp, &
		& 2.57099E+00_rp, 0.221638E+00_rp, 0.235152E+00_rp, 0.256775E+00_rp, 0.302725E+00_rp, 0.397326E+00_rp, 0.589232E+00_rp, 0.973043E+00_rp, 1.74337E+00_rp, 3.29754E+00_rp, 6.42749E+00_rp]
    Group_org_j = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
    !# BAT Functional Group Options		
    !# 1: hydroxyl		
    !# 2: carbonyl		
    !# 3: ketone		
    !# 4: hydroperoxide		
    !# 5: ether		
    !# 6: ester		
    !# 7: hydroperoxideSOA		
    !# 8: SOA chemicals		
    !# 9: PEG		

    ! calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at any RH > 0.0:
    ! calculate the equilibrium water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at any RH > 0.0:
    ! calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at RH = 0.0 (these are only needed later, for the calculation of kappaHGF (see last step):
    ! calculate the equilibrium water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at RH = 0.0 (these are only needed later, for the calculation of kappaHGF (see last step):
    ! estimate the density of each organic species j (these density estimates are only needed later, for the calculation of kappaHGF (see last step)
    
    ! estimate Csat from C* at dry conditions
	CALL CSAT_APPROX(N_org, M_org_j, O2C_org_j, C_star_org_j, C_gas_particle_org_j, Group_org_j, C_sat_org_j)
		
	! run coupled BAT+VBS model at any RH >= 0.0
    CALL BAT_INVERSION_VBS(N_a_water, N_org, M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_org_j, C_gas_particle_org_j, &
        & Group_org_j, activity_water, C_PM_org_j_alpha, C_PM_org_j_beta, C_PM_water_alpha, C_PM_water_beta)	
	
	! example output to screen:
    WRITE(*,'(A,ES13.6)') "Selected predicted properties for given input at RH = ", activity_water
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in liquid phase alpha (ug/m3): ", C_PM_water_alpha + SUM(C_PM_org_j_alpha)
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in liquid phase beta (ug/m3): ", C_PM_water_beta + SUM(C_PM_org_j_beta)
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in the particle (ug/m3): ", (C_PM_water_alpha + SUM(C_PM_org_j_alpha)) + (C_PM_water_beta + SUM(C_PM_org_j_beta))
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of water in the particle (ug/m3): ", C_PM_water_alpha + C_PM_water_beta
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of organics in the particle (ug/m3): ", SUM(C_PM_org_j_alpha) + SUM(C_PM_org_j_beta)
    WRITE(*,'(A)') "Example BAT_INVERSION_VBS program completed"
    READ(*,*) ! wait for user action

END PROGRAM MAIN_PROG

