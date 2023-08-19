!****************************************************************************************
!*   PROGRAM EXAMPLE_BAT_INVERSION                                                      *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Fortran program that shows how to use the outputs of the BAT_INVERSION_SIMULATION  *
!*   subroutine.                                                                        *
!*   Equivalent to the Python example file example_bat_inversion.py.                    *
!*                                                                                      *
!*   SUROUTINE BAT_INVERSION_SIMULATION                                                 *
!*   Input : a mixture of multiple organic species (here 10 as an example) and water.   *
!*   Output:                                                                            *
!*   1) activity coefficient of each organic species j in liquid phases                 *
!*      alpha (water-rich) and beta (organic-rich),                                     *
!*   2) mass fraction of each organic species j in liquid phases                        *
!*      alpha (water-rich) and beta (organic-rich),                                     *
!*   3) total (organics + water) mass concentration (ug/m3) in liquid phases            * 
!*      alpha (water-rich) and beta (organic-rich).                                     *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 2021-04-15                                                             *
!*   -> latest changes: 2022-05-04                                                      *
!****************************************************************************************
    
PROGRAM EXAMPLE_BAT_INVERSION
    
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

    
    !On variable KIND: For compatibility with typical choices in other Fortran programs and f2py,  
    !we use INTEGER(4) and REAL(8) as precision KINDs in this main program (serving as interface). 
    !Within the BAT procedures, a set working real precision (wrp) and similar for integer is used, 
    !which are defined in module PRECISION_MOD (and accessible via BATVBS_MOD).
    
    ! Inputs
    INTEGER, PARAMETER          :: ip = KIND(1.0E0)                                                     ! Precision of INTEGERS
    INTEGER, PARAMETER          :: rp = KIND(1.0D0)                                                     ! Precision of REALS
    INTEGER(ip),PARAMETER       :: N_org = 10                                                           ! Number of organic species in the system
    INTEGER(ip),DIMENSION(N_org):: Group_org_j                                                          ! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(rp), DIMENSION(N_org)  :: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, &                          ! Molar mass (g/mol) and elemental ratios (O/C, H/C, N/C) of each organic species j in the system
                                                                                                        ! * H/C and N/C are optional inputs; set their value to a negative real number when it is not known *
                                    & C_PM_org_j, &                                                     ! "dry" mass concentration of each organic species j in the particle phase (ug/m3) (here treated as given)
                                    & C_sat_org_j                                                       ! Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
    REAL(rp)                    :: activity_water                                                       ! Target relative humidity (gas-phase water activity) that the mole fraction of water (and organics) in the mixture must match (scale is from 0 to 1)
    REAL(rp), PARAMETER         :: M_w = 18.015D0                                                       ! Molar mass of water (g/mol)
    
    ! Local variables!
    ! mass concentrations at a given RH
    REAL(rp)                     :: C_PM_water_org_alpha, C_PM_water_org_beta, &                        ! Organic + water mass concentration in liquid phases alpha and beta that results from BAT calculations (ug/m3)
                                    & C_PM_water_org_alpha_beta, &                                      ! Organic + water mass concentration in overall liquid phase alpha + beta that results from BAT calculations (ug/m3)
                                    & C_PM_org_alpha, C_PM_org_beta, &                                  ! Organic mass concentration in liquid phases alpha and beta that results from BAT calculations (ug/m3) 
                                    & C_PM_org_alpha_beta, &                                            ! Organic mass concentration in liquid phase alpha + beta that results from BAT calculations (ug/m3)                                               
                                    & C_PM_water_alpha, C_PM_water_beta                                 ! Water mass concentration in liquid phases alpha and beta that results from BAT calculations (ug/m3)
    REAL(rp), DIMENSION(N_org)   :: C_PM_org_j_alpha,  C_PM_org_j_beta, &                               ! Mass concentration of each organic species j in liquid phases alpha and beta that results from BAT calculations (ug/m3)
                                    & C_star_org_j_alpha, C_star_org_j_beta, &                          ! Instantaneous (before VBS treatment) effective saturation concentration of each organic species j via liquid phases alpha and beta that results from BAT calculations (ug/m3)
                                    & C_star_org_j                                                      ! Instantaneous (before VBS treatment) effective saturation concentration of each organic species j (weighted mean using q_org_j_alpha values as weights) that results from BAT calculations (ug/m3)
                                           
    ! activity coefficients at a given RH
    REAL(rp), DIMENSION(N_org)   :: activity_coefficient_org_j_alpha, activity_coefficient_org_j_beta   ! Activity coefficient of each organic species j in liquid phases alpha and beta that results from BAT calculations
    
    ! mass fractions
    REAL(rp)                     :: mass_fraction_water_alpha, mass_fraction_water_beta                 ! Mass fraction of water in liquid phases alpha (with respect to C_PM_water_org_alpha) and beta (with respect to C_PM_water_org_beta) that results from BAT calculations
    REAL(rp), DIMENSION(N_org)   :: mass_fraction_org_j_alpha, mass_fraction_org_j_beta                 ! Mass fraction of organic species j in liquid phases alpha (with respect to C_PM_water_org_alpha) and beta (with respect to C_PM_water_org_beta) that results from BAT calculations
    
    ! density estimates
    REAL(rp), DIMENSION(N_org)   :: densityEst_org_j_g_cm3                                              ! Density of each organic species j (g/cm3)
    REAL(rp)                     :: densityEst_water_g_cm3                                              ! Density of water (g/cm3)
    
    ! fractional liquid僕iquid partitioning at a given RH
    REAL(rp), DIMENSION(N_org)   :: q_org_j_alpha, q_org_j_beta                                         ! Fractional liquid僕iquid partitioning of each organic species j to liquid phases alpha and beta that results from BAT calculations
    
    ! volume contributions at a given RH
    REAL(rp)                     :: V_org, V_water                                           			! Cumulative contribution of organic component volumes at a given RH; cumulative contribution of water volume at a given RH (ug cm3 / m3 / g)
    
    ! kappa
    REAL(rp)                     :: kappaHGF                                                            ! Hygroscopicity parameter kappa of OA
    
    ! set the relative humidity
    activity_water = 0.40D0    ! RH = 0.40  (= 40 %)

    ! the first element corresponds to the first organic species in the system, the second element corresponds to the second organic species in the system, etc.
    M_org_j = [200.0D+00, 188.0D+00, 216.0D+00, 368.0D+00, 368.0D+00, 204.0D+00, 195.0D+00, 368.0D+00, 158.0D+00, 206.0D+00]
    O2C_org_j = [0.4D+00, 0.444D+00, 0.5D+00, 0.368D+00, 0.368D+00, 0.556D+00, 0.857D+00, 0.368D+00, 0.375D+00, 0.75D+00 ]  ! O/C is not optional
    H2C_org_j = [1.6D+00, 1.78D+00, 1.6D+00, 1.47D+00, 1.56D+00, 1.78D+00, 1.75D+00, 1.56D+00, 1.75D+00, 1.75D+00]          ! H/C is optional; set H/C to a negative real number if it is not known
    N2C_org_j = [-1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00]        ! N/C is optional; set N/C to a negative real number if it is not known 
    C_sat_org_j = [5.74D+03, 3.27D+02, 1.67D+02, 2.79D-06, 1.05D+02, 2.13D+00, 7.19D-01, 3.64D-06, 1.16D+03, 3.02D-02]
    C_PM_org_j  = [4.79D+00 , 1.98D+00 , 0.86D+00 , 2.07D+00 , 0.72D+00 , 0.519D+00, 0.366D+00, 0.82D+00, 0.199D+00, 0.213D+00]
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
    
    ! calculate the activity coefficient of each organic species j in liquid phases alpha and beta that results from BAT calculations at a given RH, 
    ! calculate the mass fraction of organic species j in liquid phases alpha (with respect to C_PM_water_org_alpha) and beta (with respect to C_PM_water_org_beta) that results from BAT calculations at a given RH,
    ! calculate the organic + water mass concentration in liquid phases alpha and beta that results from BAT calculations at a given RH:   
    ! estimate the density of each organic species j (these density estimates are only needed later, for the calculation of kappaHGF (see last step) 
    BLOCK
        USE BATVBS_MOD ! Module that contains the main BAT and VBS subroutines
        
        REAL(wrp)                   :: activity_water_wrp, C_PM_water_org_alpha_wrp, C_PM_water_org_beta_wrp
        REAL(wrp), DIMENSION(N_org) :: M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, C_sat_org_j_wrp, C_PM_org_j_wrp, &
                                        & activity_coefficient_org_j_alpha_wrp, activity_coefficient_org_j_beta_wrp, &
                                        & mass_fraction_org_j_alpha_wrp, mass_fraction_org_j_beta_wrp, densityEst_org_j_g_cm3_wrp

        ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
        activity_water_wrp = REAL(activity_water, kind = wrp)
        M_org_j_wrp = REAL(M_org_j, kind = wrp)
        O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
        H2C_org_j_wrp = REAL(H2C_org_j, kind = wrp)
        N2C_org_j_wrp = REAL(N2C_org_j, kind = wrp)
        C_sat_org_j_wrp = REAL(C_sat_org_j, kind = wrp)
        C_PM_org_j_wrp = REAL(C_PM_org_j, kind = wrp)
        densityEst_org_j_g_cm3_wrp = REAL(densityEst_org_j_g_cm3, kind = wrp)
        activity_coefficient_org_j_alpha_wrp = REAL(activity_coefficient_org_j_alpha, kind = wrp)
        activity_coefficient_org_j_beta_wrp = REAL(activity_coefficient_org_j_beta, kind = wrp)
        mass_fraction_org_j_alpha_wrp = REAL(mass_fraction_org_j_alpha, kind = wrp)
        mass_fraction_org_j_beta_wrp = REAL(mass_fraction_org_j_beta, kind = wrp)
        C_PM_water_org_alpha_wrp = REAL(C_PM_water_org_alpha, kind = wrp)
        C_PM_water_org_beta_wrp = REAL(C_PM_water_org_beta, kind = wrp)
		
        CALL BAT_INVERSION_SIMULATION(N_org, M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, Group_org_j, C_PM_org_j_wrp, &
            & activity_water_wrp, activity_coefficient_org_j_alpha_wrp, activity_coefficient_org_j_beta_wrp, mass_fraction_org_j_alpha_wrp, &
            & mass_fraction_org_j_beta_wrp, C_PM_water_org_alpha_wrp, C_PM_water_org_beta_wrp)    
        
        ! Org_density_Estimate_KGv subroutine estimates the density of each organic species j; these density estimates are only needed later, for the calculation of kappaHGF (see last step):
        CALL Org_density_Estimate_KGv1(M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, densityEst_org_j_g_cm3_wrp)
    
        ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
        activity_coefficient_org_j_alpha = REAL(activity_coefficient_org_j_alpha_wrp, kind = rp)
        activity_coefficient_org_j_beta = REAL(activity_coefficient_org_j_beta_wrp, kind = rp)
        mass_fraction_org_j_alpha = REAL(mass_fraction_org_j_alpha_wrp, kind = rp)
        mass_fraction_org_j_beta = REAL(mass_fraction_org_j_beta_wrp, kind = rp)
        C_PM_water_org_alpha = REAL(C_PM_water_org_alpha_wrp, kind = rp)
        C_PM_water_org_beta = REAL(C_PM_water_org_beta_wrp, kind = rp)
        densityEst_org_j_g_cm3 = REAL(densityEst_org_j_g_cm3_wrp, kind = rp)
    END BLOCK
    
    ! calculate the mass concentration of each organic species j in liquid phases alpha and beta that results from BAT calculations at a given RH:
    C_PM_org_j_alpha = C_PM_water_org_alpha * mass_fraction_org_j_alpha
    C_PM_org_j_beta = C_PM_water_org_beta * mass_fraction_org_j_beta
    
    ! calculate the (total) organic mass concentration in liquid phases alpha and beta that results from BAT calculations at a given RH:
    C_PM_org_alpha = SUM(C_PM_org_j_alpha)
    C_PM_org_beta = SUM(C_PM_org_j_beta) 
    
    ! calculate the mass fraction of water in liquid phases alpha (with respect to C_PM_water_org_alpha) and 
    !beta (with respect to C_PM_water_org_beta) that results from BAT calculations at a given RH:
    mass_fraction_water_alpha = 1.0D+00 - SUM(mass_fraction_org_j_alpha)
    mass_fraction_water_beta = 1.0D+00 - SUM(mass_fraction_org_j_beta)
        
    ! calculate the water mass concentration in liquid phases alpha and beta that results from BAT calculations at a given RH:
    C_PM_water_alpha = C_PM_water_org_alpha * mass_fraction_water_alpha
    C_PM_water_beta = C_PM_water_org_beta * mass_fraction_water_beta
    
    ! calculate the organic + water mass concentration in liquid phase alpha + beta that results from BAT calculations at a given RH:
    C_PM_water_org_alpha_beta = C_PM_water_org_alpha + C_PM_water_org_beta
    ! or C_PM_water_org_alpha_beta = (C_PM_org_alpha + C_PM_water_alpha) + (C_PM_org_beta + C_PM_water_beta)
    
    ! calculate the fractional liquid僕iquid partitioning of each organic species j to liquid phases alpha and 
    !beta that results from BAT calculations at a given RH:
    ! (in the case of a liquid僕iquid equilibrium, the relative phase preferences are described by q_org_j_alpha and q_org_j_beta)
    q_org_j_alpha = C_PM_org_j_alpha/(C_PM_org_j_alpha + C_PM_org_j_beta)  ! or C_PM_org_j_alpha/(C_PM_org_j)
    q_org_j_beta = 1.0D0 - q_org_j_alpha
    
    ! calculate the instantaneous (before VBS treatment) effective saturation concentration C* according to equation 6 in https://doi.org/10.5194/acp-19-13383-2019:
    ! via liquid phase alpha
    C_star_org_j_alpha = C_sat_org_j * activity_coefficient_org_j_alpha * q_org_j_alpha * C_PM_water_org_alpha_beta/(M_org_j * (SUM(C_PM_org_j_alpha/(M_org_j)) + (C_PM_water_alpha/M_w)))
    ! via liquid phase beta
    C_star_org_j_beta = C_sat_org_j * activity_coefficient_org_j_beta * q_org_j_beta * C_PM_water_org_alpha_beta/(M_org_j * (SUM(C_PM_org_j_beta/(M_org_j)) + (C_PM_water_beta/M_w)))    
    ! weighted mean using q_org_j_alpha values as weights
    C_star_org_j = (C_star_org_j_alpha * q_org_j_alpha) + (C_star_org_j_beta * q_org_j_beta)

    ! estimate the hygroscopicity parameter kappa of the OA:
    ! this is an estimation of kappa of the OA based on the water uptake predicted by the BAT model (i.e., before VBS).
    
    ! a) calculate the volume contributions:
    densityEst_water_g_cm3 = 0.997D0                                                      ! density of water in g/cm3
    
    V_org = SUM((C_PM_org_j_alpha + C_PM_org_j_beta)/densityEst_org_j_g_cm3)              ! cumulative contribution of organic component volumes at a given RH in ug cm3 / m3 / g
    V_water = (C_PM_water_alpha + C_PM_water_beta)/densityEst_water_g_cm3                 ! cumulative contribution of water volume at a given RH in ug cm3 / m3 / g
    
    ! b) calculate the hygroscopicity parameter kappa of OA:
    kappaHGF = ((1.0D0/activity_water)-1.0D0)*(V_water/V_org)
    
    ! example output to screen:
    WRITE(*,'(A,ES13.6)') "Selected predicted properties for given input at RH = ", activity_water
    WRITE(*,'(A,ES13.6)') "Mass concentration of org and water in liquid phase alpha (ug/m3): ", C_PM_water_org_alpha
    WRITE(*,'(A,ES13.6)') "Mass concentration of org and water in liquid phase beta (ug/m3): ", C_PM_water_org_beta
    WRITE(*,'(A,ES13.6)') "Mass concentration of org and water in the particle (phase alpha + beta) (ug/m3): ", C_PM_water_org_alpha_beta    
    WRITE(*,'(A,*(ES13.6,",",1X))') "Instantaneous effective saturation concentration, C*_j, of each organic species j (ug/m3): ", C_star_org_j
    WRITE(*,'(A,ES13.6)') "Hygroscopicity parameter kappa of the organic aerosol: ", kappaHGF 
    WRITE(*,'(A)') "Example BAT_INVERSION program completed"
    READ(*,*) ! wait for user action
    
END PROGRAM EXAMPLE_BAT_INVERSION