!****************************************************************************************
!*   PROGRAM EXAMPLE_BAT_INVERSION_VBS                                                  *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Fortran program that shows how to use the outputs of the                           *
!*   BAT_INVERSION_VBS_SIMULATION subroutine.                                           *
!*   Equivalent to the Python example file example_bat_inversion_vbs.py.                *
!*                                                                                      *
!*   SUROUTINE BAT_INVERSION_VBS_SIMULATION                                             *
!*   Input system: A mixture of multiple organic species and water.                     *
!*   Output:                                                                            *
!*   1) equilibrium mass concentration (ug/m3) of each organic species j in liquid      *
!*      phases alpha (water-rich) and beta (organic-rich) that results from BAT+VBS     *
!*      calculations (ug/m3),                                                           *
!*   2) equilibrium water mass concentration (ug/m3) in liquid phases alpha (water-rich)* 
!*      and beta (organic-rich) that results from BAT+VBS calculations.                 *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 2021-04-15                                                             *
!*   -> latest changes: 2022-05-04                                                       *
!****************************************************************************************
    
PROGRAM EXAMPLE_BAT_INVERSION_VBS
        
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

    
    ! Inputs
    INTEGER, PARAMETER                      :: ip = KIND(1.0E0)                                         ! Precision of INTEGERS
    INTEGER, PARAMETER                      :: rp = KIND(1.0D0)                                         ! Precision of REALS
    INTEGER(ip),PARAMETER                   :: N_org = 10                                               ! Number of organic species in the system
    INTEGER(ip), DIMENSION(N_org)           :: Group_org_j                                              ! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(rp), DIMENSION(N_org)              :: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, &              ! Molecular weight (g/mol) and elemental ratios (O/C, H/C, N/C) of each organic species j in the system
                                                                                                        ! * H/C and N/C are optional inputs; set their value to a negative real number when it is not known *
                                                & C_gas_particle_org_j, &                               ! total mass concentration of each organic species j in gas phase + particle phase (ug/m3)
                                                & C_sat_org_j                                           ! Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
    REAL(rp)                                :: activity_water, activity_water_dry                       ! Target relative humidity (gas-phase water activity) that the mole fraction of water (and organics) in the mixture must match (scale is from 0 to 1)

    ! Local variables
    ! mass concentrations at a given RH
    REAL(rp)                                :: C_PM_water_alpha, C_PM_water_beta, &                     ! Equilibrium water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH (ug/m3)
                                                & C_PM_water_org_alpha, C_PM_water_org_beta, &          ! Equilibrium organic + water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH (ug/m3) 
                                                & C_PM_water_org_alpha_beta, &                          ! Equilibrium organic + water mass concentration in liquid phase alpha + beta that results from BAT+VBS calculations at a given RH (ug/m3)
                                                & C_PM_org_alpha_beta                                   ! Equilibrium organic mass concentration in liquid phase alpha + beta that results from BAT+VBS calculations at a given RH (ug/m3)
    REAL(rp), DIMENSION(N_org)              :: C_PM_org_j_alpha, C_PM_org_j_beta                        ! Equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH (ug/m3)
    
    ! density estimates
    REAL(rp), DIMENSION(N_org)              :: densityEst_org_j_g_cm3                                   ! Density of each organic species j (g/cm3)
    REAL(rp)                                :: densityEst_water_g_cm3                                   ! Density of water (g/m3)

    ! volume contributions at a given RH
    REAL(rp)                     			:: V_org, V_water                                           ! Cumulative contribution of organic component volumes at a given RH; cumulative contribution of water volume at a given RH (ug cm3 / m3 / g)
    
    ! kappa
    REAL(rp)                                :: kappaHGF                                                 ! Hygroscopicity parameter kappa of OA, related to the hygroscopic growth factor of the organic mixture as a function of composition (and indirectly RH)
    
    ! set relative humidity values
    activity_water = 0.40D0         ! equivalent to RH = 0.85 = 85 %
     
    ! the first element corresponds to the first organic species in the system, the second element corresponds to the second organic species in the system, etc
    M_org_j = [200.0D+00, 188.0D+00, 216.0D+00, 368.0D+00, 368.0D+00, 204.0D+00, 195.0D+00, 368.0D+00, 158.0D+00, 206.0D+00]
    O2C_org_j = [0.4D+00, 0.444D+00, 0.5D+00, 0.368D+00, 0.368D+00, 0.556D+00, 0.857D+00, 0.368D+00, 0.375D+00, 0.75D+00 ]  ! O/C is not optional
    H2C_org_j = [1.6D+00, 1.78D+00, 1.6D+00, 1.47D+00, 1.56D+00, 1.78D+00, 1.75D+00, 1.56D+00, 1.75D+00, 1.75D+00]          ! H/C is optional; set H/C to a negative real number if it is not known
    N2C_org_j = [-1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00]        ! N/C is optional; set N/C to a negative real number if it is not known 
    C_sat_org_j = [5.74D+03, 3.27D+02, 1.67D+02, 2.79D-06, 1.05D+02, 2.13D+00, 7.19D-01, 3.64D-06, 1.16D+03, 3.02D-02]
    C_gas_particle_org_j  = [8.79D+00 , 3.98D+00 , 1.13D+00 , 4.07D+00 , 1.02D+00 , 0.919D+00, 0.766D+00, 1.02D+00 , 0.399D+00, 0.313D+00]
    Group_org_j = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    !# BAT Functional Group Options		
    !# 1: hydroxyl		
    !# 2: carboxyl		
    !# 3: ketone		
    !# 4: hydroperoxide		
    !# 5: ether		
    !# 6: ester		
    !# 7: hydroperoxideSOA		
    !# 8: SOA chemicals		
    !# 9: PEG		
    
    ! calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH:
    ! calculate the equilibrium water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH:
    ! calculate the equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH (these are only needed later, for the calculation of kappaHGF (see last step)):
    ! estimate the density of each organic species j (these density estimates are only needed later, for the calculation of kappaHGF (see last step)
    
    BLOCK
        USE BATVBS_MOD ! Module that contains the main BAT and VBS subroutines
        
        REAL(wrp)                       ::  activity_water_wrp, C_PM_water_alpha_wrp, C_PM_water_beta_wrp
        REAL(wrp), DIMENSION(N_org)     ::  M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, C_sat_org_j_wrp, C_gas_particle_org_j_wrp, &
                                                & C_PM_org_j_alpha_wrp, C_PM_org_j_beta_wrp, densityEst_org_j_g_cm3_wrp
        REAL(wrp)                       ::  activity_water_dry_wrp, C_PM_water_alpha_dry_wrp, C_PM_water_beta_dry_wrp
        REAL(wrp), DIMENSION(N_org)     ::  C_PM_org_j_alpha_dry_wrp, C_PM_org_j_beta_dry_wrp

        ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
        activity_water_wrp = REAL(activity_water, kind = wrp)
        activity_water_dry_wrp = REAL(activity_water_dry, kind = wrp)
        M_org_j_wrp = REAL(M_org_j, kind = wrp)
        O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
        H2C_org_j_wrp = REAL(H2C_org_j, kind = wrp)
        N2C_org_j_wrp = REAL(N2C_org_j, kind = wrp)
        C_sat_org_j_wrp = REAL(C_sat_org_j, kind = wrp)
        C_gas_particle_org_j_wrp = REAL(C_gas_particle_org_j, kind = wrp)
        densityEst_org_j_g_cm3_wrp = REAL(densityEst_org_j_g_cm3, kind = wrp)
        C_PM_org_j_alpha_wrp = REAL(C_PM_org_j_alpha, kind = wrp)
        C_PM_org_j_beta_wrp = REAL(C_PM_org_j_beta, kind = wrp)
        C_PM_water_alpha_wrp = REAL(C_PM_water_alpha, kind = wrp)
        C_PM_water_beta_wrp = REAL(C_PM_water_beta, kind = wrp)

        ! run coupled BAT+VBS model at a given RH:
        CALL BAT_INVERSION_VBS_SIMULATION(N_org, M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, C_sat_org_j_wrp, C_gas_particle_org_j_wrp, &
            & Group_org_j, activity_water_wrp, C_PM_org_j_alpha_wrp, C_PM_org_j_beta_wrp, C_PM_water_alpha_wrp, C_PM_water_beta_wrp)
        
        ! Org_density_Estimate_KGv subroutine estimates the density of each organic species j; these density estimates are only needed later, for the calculation of kappaHGF (see last step):  
        CALL Org_density_Estimate_KGv1(M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, densityEst_org_j_g_cm3_wrp)
        
        ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
        C_PM_org_j_alpha = REAL(C_PM_org_j_alpha_wrp, kind = rp)
        C_PM_org_j_beta = REAL(C_PM_org_j_beta_wrp, kind = rp)
        C_PM_water_alpha = REAL(C_PM_water_alpha_wrp, kind = rp)
        C_PM_water_beta = REAL(C_PM_water_beta_wrp, kind = rp)
        densityEst_org_j_g_cm3 = REAL(densityEst_org_j_g_cm3_wrp, kind = rp)
    END BLOCK
  
    ! calculate the equilibrium organic + water mass concentration in liquid phases alpha and beta that results from BAT+VBS calculations at a given RH:
    C_PM_water_org_alpha = SUM(C_PM_org_j_alpha) + C_PM_water_alpha
    C_PM_water_org_beta = SUM(C_PM_org_j_beta) + C_PM_water_beta
    
    ! calcualte the equilibrium organic + water mass concentration in liquid phase alpha + beta that results from BAT+VBS calculations at a given RH:
    C_PM_water_org_alpha_beta = C_PM_water_org_alpha + C_PM_water_org_beta
    
    ! predict the hygroscopicity parameter kappa of OA(this is an estimation of kappa of the OA based on the outputs of the coupled BAT+VBS model):
    
    ! a) calculate the volume contributions:
    densityEst_water_g_cm3 = 0.997D+00                                                      ! density of water in g/cm3
    V_org =   SUM((C_PM_org_j_alpha + C_PM_org_j_beta)/densityEst_org_j_g_cm3)              ! cumulative contribution of organic component volumes at a given RH
    V_water = (C_PM_water_alpha + C_PM_water_beta)/densityEst_water_g_cm3                   ! cumulative contribution of water volume at a given RH
    
    ! b) calculate the hygroscopicity parameter kappa of OA:
    kappaHGF = ((1.0D0/activity_water)-1.0D0)*(V_water/V_org)
    
    ! example output to screen:
    WRITE(*,'(A,ES13.6)') "Selected predicted properties for given input at RH = ", activity_water
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in liquid phase alpha (ug/m3): ", C_PM_water_alpha + SUM(C_PM_org_j_alpha)
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in liquid phase beta (ug/m3): ", C_PM_water_beta + SUM(C_PM_org_j_beta)
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of org + water in the particle (ug/m3): ", C_PM_water_org_alpha_beta
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of water in the particle (ug/m3): ", C_PM_water_alpha + C_PM_water_beta
    WRITE(*,'(A,ES13.6)') "Equilibrium mass concentration of organics in the particle (ug/m3): ", SUM(C_PM_org_j_alpha) + SUM(C_PM_org_j_beta)
    WRITE(*,'(A,ES13.6)') "Hygroscopicity parameter kappa of the organic aerosol: ", kappaHGF 
    WRITE(*,'(A)') "Example BAT_INVERSION_VBS program completed"
    READ(*,*) ! wait for user action

END PROGRAM EXAMPLE_BAT_INVERSION_VBS