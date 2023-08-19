!****************************************************************************************
!*   MODULE BATVBS_MOD                                                                  *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that contains all the subroutines related to the BAT and BATVBS simulations.*
!*                                                                                      *
!*   This module includes the BAT_inversion_simulation subroutine, which estimates the  *
!*   (mole-fraction-based) activity coefficient of organic species at a given RH when   *
!*   the mole fraction of organics and water in each binary mixture are not known.      *
!*                                                                                      *
!*   This module includes the BAT_simulation subroutine, which estimates the            *
!*   (mole-fraction-based) activity coefficient of organic species at a given RH when   *
!*   the mole fraction of organics and water in each binary mixture are known.          *
!*                                                                                      *
!*   This module includes the BAT_INVERSION_VBS_simulation subroutine, which represents * 
!*   the complete coupled BAT+VBS evaluation that estimates the equilibrium organic     *
!*   aerosol mass concentration.                                                        *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 20-10-2022                                                      *
!****************************************************************************************

MODULE BATVBS_MOD

    USE DATATYPE_MOD                        ! Module that defines derived data types
    USE TOOLS_MOD                           ! Module that contains subroutines to perform different recurrent calculations
    USE PRECISION_MOD                       ! Module that defines the precision/length of REAL/INTERGER/CHAR data type
    USE NN_MOD                              ! Module that contains subroutines related to the neural network simulations
    USE BAT_MOD                             ! Module that contains all the constants used in the neural network simulatioLLEns
    USE MINPACK_MOD                         ! Module that contains the VBS solver
    USE OBJECTIVE_FUNCTION_VARIABLES_MOD    ! Module that contains objective function variables that are updated at each interation
    USE OBJECTIVE_FUNCTION_MOD              ! Module that contains the VBS solver objective function

    IMPLICIT NONE

CONTAINS

    !****************************************************************************************
    !*   SUBROUTINE BAT_SIMULATION                                                          *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that calculates the (mole-fraction-based) activity coefficient of       *
    !*   water and one (1) organic species in their binary mixture.                         *
    !*   (Gorkowski, K., Preston, T. C., and Zuend, A.: Relative-humidity-dependent organic *
    !*   aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem.     *
    !*   Phys., 19, 13383–13407, https://doi.org/10.5194/acp-19-13383-2019, 2019.).         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano Damha                                                               *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************
    
    SUBROUTINE BAT_SIMULATION(M_org, O2C_org, H2C_org, N2C_org, Group_org, x_org, y_water, y_org)
    
        IMPLICIT NONE
    
        ! Input
        INTEGER(wip), INTENT(IN)            :: Group_org                    ! Type of oxygen-bearing functional group of the organic species in the binary mixture
        REAL(wrp), INTENT(IN)               :: O2C_org, H2C_org, N2C_org, & ! Elemental ratios (O/C, H/C, N/C) of the organic species in the binary mixture
                                               & M_org, x_org               ! Molecular weight (g/mol) and mole fraction of the organic species in the binary mixture
        
        ! Output
        REAL(wrp), INTENT(OUT)              :: y_water, y_org               ! Activity cofficient of water and organic species in the binary mixture
                
        ! Local
        REAL(wrp), DIMENSION(1)             :: O2C_value, H2C_value, N2C_value_densityOnly, Molar_mass_ratio, &
                                               & func1, func2, Gibbs_RT, dGibbs_RTdx2, mole_frac_org, a_water, a_org, w_water, w_org, &
                                               & y_org_out, y_water_out
        CHARACTER(clen), DIMENSION(1)       :: BAT_functional_group
    
        ! Extracting properties
        O2C_value = O2C_org                 ! O/C ratio of the organic species in the binary system
        H2C_value = H2C_org                 ! H/C ratio of the organic species in the binary system
        N2C_value_densityOnly = N2C_org     ! N/C ratio of the organic species in the binary system
        mole_frac_org = x_org               ! Mole fraction of the organic species in the binary system
        Molar_mass_ratio = [M_water/M_org]  ! Molecular weight ratio M_water/M_org

        ! Set a value (== 0.0_wrp) for N/C ratio if it is not known at input
        WHERE (N2C_value_densityOnly < -tinynumber) N2C_value_densityOnly = 0.0_wrp
        
        ! Oxygen-bearing functional group of the organic species in the binary system
        IF (Group_org == 1) THEN
            BAT_functional_group(1) = 'hydroxyl'
        ELSE IF (Group_org == 2) THEN
            BAT_functional_group(1) = 'carboxyl'
        ELSE IF (Group_org == 3) THEN
            BAT_functional_group(1) = 'ketone'
        ELSE IF (Group_org == 4) THEN
            BAT_functional_group(1) = 'hydroperoxide'
        ELSE IF (Group_org == 5) THEN
            BAT_functional_group(1) = 'ether'
        ELSE IF (Group_org == 6) THEN
            BAT_functional_group(1) = 'ester'
        ELSE IF (Group_org == 7) THEN
            BAT_functional_group(1) = 'hydroperoxideSOA'
        ELSE IF (Group_org == 8) THEN
            BAT_functional_group(1) = 'SOA chemicals'
        ELSE IF (Group_org == 9) THEN
            BAT_functional_group(1) = 'PEG'       
        ELSE
            STOP
        END IF
        
        ! Calculate the activity coefficient of the organic species (1) and water in the binary mixture
        CALL BAT_properties_calculation_v1(mole_frac_org, O2C_value, H2C_value, Molar_mass_ratio, BAT_functional_group, &
            & N2C_value_densityOnly, func1, func2, y_water_out, y_org_out, a_water, a_org, w_water, &
            & w_org, Gibbs_RT, dGibbs_RTdx2) 
        
        y_water = y_water_out(1)
        y_org   = y_org_out(1)
          
    END SUBROUTINE BAT_SIMULATION
                    
    !****************************************************************************************
    !*   SUBROUTINE BAT_INVERSION_SIMULATION                                                *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that calculates (mole-fraction-based) activity coefficient of organics, *
    !*   fractional liquid-liquid partitioning of organics, mass fraction of organics,      *
    !*   mass fraction of water, mole fraction of organics and mole fraction of water       *
    !*   in phases alpha (water-rich) and beta (water-poor) at a given water activity (RH)  *
    !*   using the BAT (INVERSION) model                                                    *
    !*   (Gorkowski, K., Preston, T. C., and Zuend, A.: Relative-humidity-dependent organic *
    !*   aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem.     *
    !*   Phys., 19, 13383–13407, https://doi.org/10.5194/acp-19-13383-2019, 2019.).         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano Damha                                                               *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************
    
    SUBROUTINE BAT_INVERSION_SIMULATION(n_species, M_j, O2C_j, H2C_j, N2C_j, Group_j, C_PM_j, a_water, &
        & actcoeff_j_alpha, actcoeff_j_beta, w_j_alpha, w_j_beta, C_PM_alpha, C_PM_beta)
    
        IMPLICIT NONE
    
        ! Input
        INTEGER(wip), INTENT(IN)                       :: n_species                                 ! Number of organic species in the system
        INTEGER(wip), DIMENSION(n_species), INTENT(IN) :: Group_j                                   ! Type of oxygen-bearing functional group of each org. species
        REAL(wrp), INTENT(IN)                          :: a_water                                   ! Water activity (Relative humidity)
        REAL(wrp), DIMENSION(n_species), INTENT(IN)    :: M_j, &                                    ! Molecular weight of each org. species in the system
                                                          & O2C_j, H2C_j, N2C_j, &                  ! Elemental ratios (O/C, H/C, N/C) of each org. species in the system
                                                          & C_PM_j                                  ! Total "dry" mass concentration of each organic species in the particle (ug/m3) 
        ! Output
        REAL(wrp), INTENT(OUT)                         :: C_PM_alpha, C_PM_beta                     ! Total (all organics + water) mass concentration in liquid phases alpha and beta (ug/m3)
        REAL(wrp), DIMENSION(n_species), INTENT(OUT)   :: actcoeff_j_alpha, actcoeff_j_beta, &      ! Activity coefficient of org. in phases alpha and beta
                                                          & w_j_alpha, w_j_beta                     ! Mass fraction of organic species in phases alpha and beta
        ! Local variables
        TYPE(input_data)                      :: simulation
        TYPE(growth_input)                    :: growth
        TYPE(VBSBAT_options_input)            :: VBSBAT_options
        REAL(wrp)                             :: Caq_onlyalpha, Caq_onlybeta, Coa_onlyalpha, Coa_onlybeta, Coaaq_onlyalpha, Coaaq_onlybeta, &
                                                 & massweighted_molar_weight_alpha, massweighted_molar_weight_beta
        REAL(wrp), DIMENSION(n_species)       :: N2C_values_densityOnly, Molar_mass_ratios, O2C_eqv, molarmass_ratio_eqv, &
                                                 & O2C_single_phase_cross_point_a, O2C_single_phase_cross_point, mean_prop_mask, &
                                                 & a_w_sep_point, aw_vec, mole_frac_org_alpha, mole_frac_org_beta, mass_fraction_org_alpha, &
                                                 & mass_fraction_org_beta, mass_fraction_water_alpha, ycalc_org_alpha, ycalc_org_beta, &
                                                 & mass_fraction_water_beta, guess_partition_coefficients, partition_coefficients_temp, &
                                                 & partition_coefficients, Cstar_j_temp, mass_inPM_temp, &
                                                 & mass_fraction_inPM, q_alpha_molefrac_phase_split_org, Mratio_temp, O2C_temp, &
                                                 & q_alpha_values, Csat_j_value, C_OM_ugPm3, O2C_values, H2C_values, Molecular_weight, &
                                                 & Coa_j_onlyalpha, Coa_j_onlybeta, q_org_alpha, w_water_alpha, w_water_beta, &
                                                 & Cstar_j_via_onlyalpha, Cstar_j_via_onlybeta
        REAL(wrp), DIMENSION(1)               :: aw_series, C_OA_out, Coa_alpha, Coa_beta, q_alpha_water_out, weight_q_alpha, &
                                                 & a_w_sep_point_of_meanPM, ycalc_org_beta_temp, mass_fraction_water_beta_temp, &
                                                 & ycalc_org_alpha_temp, mass_fraction_water_alpha_temp, fit_exit_flag_save, error_save, &
                                                 & C_OA_ratio,kappa, mole_a, mole_b, O2C, H2C, Mratio, N2C, func1, func2, ycal_water, ycalc_org_a, &
                                                 & ycalc_org_b, activity_water, activity_calc2_alpha, activity_calc2_beta, &
                                                 & mass_fraction_water_a, mass_fraction_water_b, mass_fraction_org_a, mass_fraction_org_b, Gibbs_RT, &
                                                 & dGibbs_RTdx2, mole_frac_fit_a, mole_frac_fit_b, error_out, guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, aw, C_OA_PM, &
                                                 & fit_exit_flag, C_OA_temp, q_alpha_water_temp, q_alpha_water, guess_C_OA_ugPm3, &
                                                 & O2C_temp_2, Mratio_temp_2, H2C_temp_2, aw_point_mean, weight_q_a, mole_frac_org_a, &
                                                 & mole_frac_org_b
        REAL(wrp), DIMENSION(1,n_species)     :: a_w_sep_matrix, aw_series_matrix, q_alpha_vsRH_values, max_q_alpha, min_q_alpha, &
                                                 & not_min_q_alpha, partition_coefficients_out, Cstar_j_out, Coa_j_alpha, Coa_j_beta, &
                                                 & Caq_j_alpha, Caq_j_beta, gamma_alpha, gamma_beta, mass_fraction_water_alpha_out, &
                                                 & mass_fraction_water_beta_out, Coa_j_PM, Caq_j_PM
        REAL(wrp), DIMENSION(n_species+1)     :: Molecular_weight_withwater
        REAL(wrp), DIMENSION(2)               :: mole_frac_bounds_alpha, mole_frac_bounds_beta, Coa_AB, Caq_AB, Coa_AB_temp, Caq_AB_temp
        REAL(wrp), DIMENSION(2,n_species)     :: activity_coefficient_AB, mass_fraction_water_AB, Coa_j_AB, Caq_j_AB, &
                                                 & Coa_j_AB_temp, Caq_j_AB_temp
        REAL(wrp), DIMENSION(1,1)             :: weight_1_1, aw_1_1, aw_point_1_1
        INTEGER(wip)                          :: s_i, i
        CHARACTER(clen), DIMENSION(1)         :: BAT_group, BAT_refinement_mode_temp, sim_name, BAT_refinement_mode
        CHARACTER(clen), DIMENSION(n_species) :: BAT_functional_group
    
        ! Extracting properties
        O2C_values = O2C_j                  ! O/C ratio of every organic species in the system
        H2C_values = H2C_j                  ! H/C ratio of every organic species in the system
        N2C_values_densityOnly = N2C_j      ! N/C ratio of every organic species in the system
        Molecular_weight = M_j              ! Molecular weight of every organic species in the system

        ! Set a value (== 0.0_wrp) for N/C ratio if it is not known at input
        WHERE (N2C_values_densityOnly < -tinynumber) N2C_values_densityOnly = 0.0_wrp
        
        ! Oxygen-bearing functional group of every organic species in the system
        DO i = 1, n_species
            IF (Group_j(i) == 1) THEN
                BAT_functional_group(i) = 'hydroxyl'
            ELSE IF (Group_j(i) == 2) THEN
                BAT_functional_group(i) = 'carboxyl'
            ELSE IF (Group_j(i) == 3) THEN
                BAT_functional_group(i) = 'ketone'
            ELSE IF (Group_j(i) == 4) THEN
                BAT_functional_group(i) = 'hydroperoxide'
            ELSE IF (Group_j(i) == 5) THEN
                BAT_functional_group(i) = 'ether'
            ELSE IF (Group_j(i) == 6) THEN
                BAT_functional_group(i) = 'ester'
            ELSE IF (Group_j(i) == 7) THEN
                BAT_functional_group(i) = 'hydroperoxideSOA'
            ELSE IF (Group_j(i) == 8) THEN
                BAT_functional_group(i) = 'SOA chemicals'
            ELSE IF (Group_j(i) == 9) THEN
                BAT_functional_group(i) = 'PEG'       
            ELSE
                STOP
            END IF
        END DO
        
        aw_series(1) = a_water ! Water activity to run the model at
        
        ! Set simulation options (see SUBROUTINE VBSBAT_setoptions)
        CALL VBSBAT_setoptions(simulation)
        BAT_refinement_mode = simulation%BAT_refinement_mode ! BAT refinement mode
        VBSBAT_options = simulation%VBSBAT_options ! BAT and BAT+VBS simulation options/settings
    
        Molar_mass_ratios = M_water/Molecular_weight ! Molecular weight ratios
    
        ! Main step: check for phase separation conditions and calculation of q_alpha values (fractional liquid-liquid partitioning
        !         of components to phase alpha)
    
        ! Calculate new O/C and Molar mass values using shift fit data to convert functionalized molecules to a hypothetical
        ! OH-equivalent molecule
        CALL convert_chemical_structure_to_OH_eqv_v3(O2C_values, Molar_mass_ratios, BAT_functional_group, O2C_temp, &
            & Mratio_temp)
    
        ! Calculate O/C single phase cross point (O/C at the miscibility limit line as a function of molar mass)
        CALL single_phase_O2C_point_KGv3(Mratio_temp, O2C_single_phase_cross_point_a)
    
        ! Mask: when the converted O/C (OH equiv.) is less than the O/C cross point for the same
        ! molar mass, we are in the less polar region => there is a separation water activity
        O2C_single_phase_cross_point = O2C_values*0.0_wrp
        WHERE (O2C_temp < O2C_single_phase_cross_point_a) O2C_single_phase_cross_point = 1.0_wrp
        mean_prop_mask = O2C_single_phase_cross_point
    
        IF (VBSBAT_options%force_phase%onePhase == 'no') THEN ! Includes two phase option for q_alpha calculation
            IF (SUM(O2C_single_phase_cross_point) > 0.0_wrp) THEN ! Has a separation water activity separation point
                                                                  ! aw_sep_point because of lower polarity
                ! Calculate the water activity separation point (a_w_sep_point) for each species
                ! The a_w_sep_point is determined using the BAT model activities and associated Gibbs energy of mixing
                ! To approximate the location and aw width over which the liquid–liquid phase separation is prescribed
                ! to occur, we ﬁrst determine a designated reference point, the so-called water activity separation point
                ! (aw_sep_point)
                CALL biphasic_to_single_phase_RH_master_v4(O2C_values, H2C_values, Molar_mass_ratios, &
                    & BAT_functional_group, a_w_sep_point)
                
                ! Isolate miscible from miscibility gap.
                a_w_sep_point = a_w_sep_point * O2C_single_phase_cross_point
                
                ! Creates a a_w_sep_point matrix to match aw_series length and chemical species width.
                a_w_sep_matrix = SPREAD(a_w_sep_point, 1, 1)
                aw_series_matrix = SPREAD(aw_series, 2, n_species)
                
                ! Calculate q_alpha as a function of water activity aw (RH) (sigmoid function to simulate transition between liquid phases)
                ! The components' fraction in phase alpha follows a smooth transition function for q_alpha with changing RH
                ! q_alpha_vsRH_values contains q_alpha values of all the species in the system as a function of RH values
                CALL q_alpha_transfer_vs_aw_calc_v1(a_w_sep_matrix, aw_series_matrix, VBSBAT_options, q_alpha_vsRH_values)
    
                ! Set threshold: start and end q_alpha values of the sigmoid function
                max_q_alpha = q_alpha_vsRH_values * 0.0_wrp
                min_q_alpha = q_alpha_vsRH_values * 0.0_wrp
                WHERE (q_alpha_vsRH_values > VBSBAT_options%q_alpha%q_alpha_bounds(1)) max_q_alpha = 1.0_wrp
                WHERE (q_alpha_vsRH_values < VBSBAT_options%q_alpha%q_alpha_bounds(2)) min_q_alpha = 1.0_wrp
    
                ! Apply q_alpha limits
                not_min_q_alpha = min_q_alpha * 0.0_wrp
                WHERE (min_q_alpha < tinynumber) not_min_q_alpha = 1.0_wrp ! where min_q_alpha == 0.0_wrp, not_min_q_alpha = 1.0_wrp
                q_alpha_vsRH_values = max_q_alpha * q_alpha_vsRH_values
                q_alpha_vsRH_values = min_q_alpha * q_alpha_vsRH_values + not_min_q_alpha
            ELSE
                q_alpha_vsRH_values = 1.0 ! Fully miscible, there is no separation water activity
            END IF
    
        ELSE IF (VBSBAT_options%force_phase%onePhase == 'beta') THEN ! Organic-rich phase only (force one phase: q_beta = 1.0_wrp)
            q_alpha_vsRH_values = 0.0_wrp
        ELSE IF (VBSBAT_options%force_phase%onePhase == 'alpha') THEN ! Water-rich phase only (force one phase: q_alpha = 1.0_wrp
            q_alpha_vsRH_values = 1.0_wrp
        ELSE
            PRINT *, 'Specify VBSBAT_options%force_phase%onePhase in SUBROUTINE VBSBAT_options'
        END IF
        
        ! matrix and vectors to save output data
        partition_coefficients_out = 0.0_wrp
        Cstar_j_out = partition_coefficients_out
        Coa_j_alpha = partition_coefficients_out
        Coa_j_beta = partition_coefficients_out
        Caq_j_alpha = partition_coefficients_out
        Caq_j_beta = partition_coefficients_out
        gamma_alpha = partition_coefficients_out
        gamma_beta = partition_coefficients_out
        mass_fraction_water_alpha_out = partition_coefficients_out
        mass_fraction_water_beta_out = partition_coefficients_out
    
        C_OA_out = 0.0_wrp
        Coa_alpha = C_OA_out
        Coa_beta = C_OA_out
        q_alpha_water_out = C_OA_out
        weight_q_alpha = C_OA_out
        a_w_sep_point_of_meanPM = C_OA_out
    
        ycalc_org_beta_temp = C_OA_out
        mass_fraction_water_beta_temp = C_OA_out
        ycalc_org_alpha_temp = C_OA_out
        mass_fraction_water_alpha_temp = C_OA_out
    
    
        s_i = 1 ! Unlike the MATLAB version of the BATVBS model, this version only accepts 1 value of RH at input
        aw(1) = aw_series(s_i) ! Select one water activity walue
        aw_vec = aw(1) ! Create a water activity vector that matches number of organic species
        
        ! Main step: use neural network (NN) to get mole fraction of organics at a given water activity (inversion: input: aw ; output: x_org)
        ! We use NN with the BAT model to ﬁnd the correct x_org input for subsequent calculations, since in most applications
        ! aw is known but not x_org
        CALL inverted_NNBAT_v8(O2C_values, H2C_values, Molar_mass_ratios, aw_vec, BAT_functional_group, &
            mole_frac_org_alpha, mole_frac_org_beta)
                
        mass_fraction_water_alpha = 0.0_wrp
        ycalc_org_alpha = mass_fraction_water_alpha
        mass_fraction_water_beta = 0.0_wrp
        ycalc_org_beta = mass_fraction_water_beta
    
        DO i = 1, n_species ! Start loop through organic species.
           
            ! Set mole fraction limits
            mole_frac_bounds_alpha = [0.0_wrp, 1.0_wrp]
            mole_frac_bounds_beta = [0.0_wrp, 1.0_wrp]
            
            ! IF (VBSBAT_options%BAT_refinement_aw >= aw(1)) THEN ! Set when the refinement should start
                ! BAT_refinement_mode_temp = ['none']  ! No interpolation is necessary below RH of BAT_refinement_aw
            ! ELSE
                ! BAT_refinement_mode_temp = BAT_refinement_mode  ! Interpolate via 501 org mole fraction points when RH is above BAT_refinement_aw
                                                                ! ! Iterative reﬁnement is required for good agreement with the targeted aw
            ! END IF
			
            BAT_refinement_mode_temp = BAT_refinement_mode
            ! Select one species of the system at a time.
            mole_frac_org_a(1) = mole_frac_org_alpha(i) ! Every x_org in phase alpha (water rich)
            mole_frac_org_b(1) = mole_frac_org_beta(i)  ! Every x_org in phase beta (organic-rich)
            O2C(1) = O2C_values(i) ! every O:C
            H2C(1) = H2C_values(i) ! every H:C
            Mratio(1) = Molar_mass_ratios(i) ! every molar mass ratio
            BAT_group(1) = BAT_functional_group(i) ! every funtional group
            N2C(1) = N2C_values_densityOnly(i) ! every N:C
            ycalc_org_a(1) = ycalc_org_alpha(i) ! every org activity coefficient in phase alpha
            ycalc_org_b(1) = ycalc_org_beta(i) ! every org activity coefficient in phase beta
            mass_fraction_water_a(1) = mass_fraction_water_alpha(i) ! every water mass fraction in phase alpha
            mass_fraction_water_b(1) = mass_fraction_water_beta(i) ! every water mass fraction in phase beta
            
            ! Use x_org to find water mass fraction and activity coefficient associated with each organic species.
            ! Iterative reﬁnement is required for good agreement with the targeted aw (VBSBAT_options%BAT_refinement_aw)
            IF (BAT_refinement_mode_temp(1) == 'interpolate') THEN ! Interpolate in both liquid phases (alpha and beta)
				
                BAT_refinement_mode_temp = ['interpolatealpha']
                CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_a, O2C, H2C, Mratio, BAT_group, &
                    & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_a, &
                    & activity_water, activity_calc2_alpha, mass_fraction_water_a, mass_fraction_org_a, Gibbs_RT, &
                    & dGibbs_RTdx2, mole_frac_fit_a, error_out)
                BAT_refinement_mode_temp = ['interpolatebeta']
                CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_b, O2C, H2C, Mratio, BAT_group, &
                    & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_b, &
                    & activity_water, activity_calc2_beta, mass_fraction_water_b, mass_fraction_org_b, Gibbs_RT, &
                    & dGibbs_RTdx2, mole_frac_fit_b, error_out)
            
			ELSE ! No refinement is necessary
				
                CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_a, O2C, H2C, Mratio, BAT_group, &
                    & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_a, &
                    & activity_water, activity_calc2_alpha, mass_fraction_water_a, mass_fraction_org_a, Gibbs_RT, &
                    & dGibbs_RTdx2, mole_frac_fit_a, error_out)
                CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_b, O2C, H2C, Mratio, BAT_group, &
                    & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_b, &
                    & activity_water, activity_calc2_beta, mass_fraction_water_b, mass_fraction_org_b, Gibbs_RT, &
                    & dGibbs_RTdx2, mole_frac_fit_b, error_out)
           
            END IF
            ! Store activity cofficients and water mass fractions (in phase alpha and beta) associated with every
            ! species in the system at a particular RH
            ycalc_org_alpha(i) = ycalc_org_a(1)
            ycalc_org_beta(i) = ycalc_org_b(1)
            mass_fraction_water_alpha(i) = 1.0_wrp - mole_frac_fit_a(1)!mass_fraction_water_a(1)
            mass_fraction_water_beta(i) = 1.0_wrp - mole_frac_fit_b(1)!mass_fraction_water_b(1)
            mass_fraction_org_alpha(i) = mole_frac_fit_a(1)!mass_fraction_org_a(1)
            mass_fraction_org_beta(i) = mole_frac_fit_a(1)!mass_fraction_org_b(1)
            
        END DO ! End loop through organic species
        
        ! y_org_alpha, y_org_beta, w_alpha, w_beta at a given aw
        activity_coefficient_AB(1,:) = ycalc_org_alpha              ! Activity coefficient (of each org) in phase alpha
        activity_coefficient_AB(2,:) = ycalc_org_beta               ! Activity coefficient (of each org) in phase beta
        mass_fraction_water_AB(1,:) = mass_fraction_water_alpha     ! Water mass fraction (per org) in phase alpha
        mass_fraction_water_AB(2,:) = mass_fraction_water_beta      ! Water mass fraction (per org) in phase beta
        
        ! Outputs
        actcoeff_j_alpha = activity_coefficient_AB(1,:)          ! Activity coeff. of organics in phases alpha and beta
        actcoeff_j_beta = activity_coefficient_AB(2,:)
        q_org_alpha = q_alpha_vsRH_values(s_i,:)    ! Fractional liquid-liquid partitioning of organics to phases alpha and beta  
        w_water_alpha = mass_fraction_water_alpha   ! Mass fraction of water in phases alpha and beta for each binary system (organic+water)
        w_water_beta = mass_fraction_water_beta
        
        ! Mass concentration calculation
        ! Phase alpha
        Coa_j_onlyalpha = C_PM_j * q_org_alpha                                          ! Organic mass in phase alpha of each organic species  
        Coa_j_onlyalpha = C_PM_j * q_org_alpha                                          ! Organic mass in phase alpha of each organic species  
        Caq_onlyalpha = SUM(Coa_j_onlyalpha*w_water_alpha/(1.0_wrp - w_water_alpha))    ! Total water mass in phase alpha based on organic mass and mass fraction of water from BAT model
        Coa_onlyalpha = SUM(Coa_j_onlyalpha)                                            ! Total organic mass in phase alpha 
        C_PM_alpha = Coa_onlyalpha+Caq_onlyalpha                                        ! Total mass (water plus organic mass) in phase alpha
    
        ! Phase beta
        Coa_j_onlybeta = C_PM_j * (1.0_wrp - q_org_alpha) ! Organic mass in phase beta of each organic species
        Caq_onlybeta = SUM(Coa_j_onlybeta * w_water_beta/(1.0_wrp - w_water_beta)) ! Total water mass in phase beta based on organic mass and mass fraction of water from BAT model
        Coa_onlybeta = SUM(Coa_j_onlybeta) ! Total organic mass in phase beta
        C_PM_beta = Coa_onlybeta + Caq_onlybeta ! Total mass (water plus organic mass) in phase beta
                
        ! Mass fractions
        w_j_alpha = Coa_j_onlyalpha/C_PM_alpha ! w_water_alpha         ! Mass fraction of organic species j that is in phase alpha
        w_j_beta = Coa_j_onlybeta/C_PM_beta ! w_water_beta             ! Mass fraction of organic species j that is in phase beta
    
    END SUBROUTINE BAT_INVERSION_SIMULATION
    
    !****************************************************************************************
    !*   SUBROUTINE BAT_INVERSION_VBS_SIMULATION                                      *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Main subroutine for the BAT+VBS simulation . There are 4 main steps:               *
    !*   1) Calculate the water activity at which liquid–liquid phase separation happens    *
    !*      for organics.                                                                   *
    !*   2) Estimate the mole fraction of organics via BAT - Neural Network.                *
    !*   3) Calculate the activity coefficient of organics at a given water activity (RH).  *
    !*   4) Use the hybrd1 equation solver to find the partitioning coefficients of organics*
    !*      and calculate the equilibrium organic aerosol mass concentration.               *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************
    
    SUBROUTINE BAT_INVERSION_VBS_SIMULATION(n_rh, n_species, M_j, O2C_j, H2C_j, N2C_j, Csat_j, Ctotal_j, Group_j, aw_bat, Coa_org_alpha, Coa_org_beta, &
        & Coa_total_water_alpha, Coa_total_water_beta)
    
        IMPLICIT NONE
    
        ! Input
		INTEGER(wip), INTENT(IN)					   :: n_rh                          ! Number of water activity/RH values
        INTEGER(wip), INTENT(IN)                       :: n_species                     ! Number of organic species in the system
        INTEGER(wip), DIMENSION(n_species), INTENT(IN) :: Group_j                       ! Type of oxygen-bearing functional group of each org. species
        REAL(wrp), DIMENSION(n_rh), INTENT(IN)		   :: aw_bat                        ! Water activity (Relative humidity)
        REAL(wrp), DIMENSION(n_species), INTENT(IN)    :: M_j, O2C_j, H2C_j, N2C_j, &   ! Molar mass, elemental ratios (O/C, H/C, N/C) of each org. species
                                                          & Csat_j, Ctotal_j            ! Total mass concentration (gas+particle) and saturation mass
                                                                                        ! concentration of each org. species       
        
        ! Output                                          
        REAL(wrp), DIMENSION(n_rh), INTENT(OUT)		         :: Coa_total_water_alpha, Coa_total_water_beta ! Total equilibrium aerosol mass concentration of water in phases alpha and beta (ug/m3)
        REAL(wrp), DIMENSION(n_rh, n_species), INTENT(OUT)   :: Coa_org_alpha, Coa_org_beta ! Equilibrium mass concentration of each organic species in phases alpha and beta (ug/m3)
        ! Local variables
        TYPE(input_data)                      :: simulation
        TYPE(specialoptions_input)            :: special_options
        TYPE(growth_input)                    :: growth
        TYPE(VBSBAT_options_input)            :: VBSBAT_options
        REAL(wrp), DIMENSION(1)               :: C_OA_PM_mat, Caq_PM_mat
        REAL(wrp), DIMENSION(n_species)       :: N2C_values_densityOnly, Molar_mass_ratios, O2C_eqv, molarmass_ratio_eqv, &
                                                 & O2C_single_phase_cross_point_a, O2C_single_phase_cross_point, mean_prop_mask, a_w_sep_point, aw_vec, &
                                                 & mole_frac_org_alpha, mole_frac_org_beta, mass_fraction_water_alpha, ycalc_org_alpha, ycalc_org_beta, &
                                                 & mass_fraction_water_beta, guess_partition_coefficients, partition_coefficients_temp, Cstar_j, &
                                                 & partition_coefficients, Cstar_j_temp, mass_inPM_temp, mass_fraction_inPM, q_alpha_molefrac_phase_split_org, &
                                                 & Mratio_temp, O2C_temp, q_alpha_values, Csat_j_value, C_OM_ugPm3, O2C_values, H2C_values, Molecular_weight
        REAL(wrp), DIMENSION(1)               :: aw_series, C_OA_out, Coa_alpha, Coa_beta, q_alpha_water_out, weight_q_alpha, &
                                                 & a_w_sep_point_of_meanPM, ycalc_org_beta_temp, mass_fraction_water_beta_temp, ycalc_org_alpha_temp, &
                                                 & mass_fraction_water_alpha_temp, fit_exit_flag_save, error_save, C_OA_ratio, kappa
        REAL(wrp), DIMENSION(1,n_species)     :: a_w_sep_matrix, aw_series_matrix, q_alpha_vsRH_values, max_q_alpha, &
                                                 & min_q_alpha, not_min_q_alpha, partition_coefficients_out, Cstar_j_out, Coa_j_alpha, Coa_j_beta, Caq_j_alpha, &
                                                 & Caq_j_beta, gamma_alpha, gamma_beta, mass_fraction_water_alpha_out, mass_fraction_water_beta_out, Coa_j_PM, &
                                                 & Caq_j_PM
        REAL(wrp), DIMENSION(n_species+1)     :: Molecular_weight_withwater
        REAL(wrp), DIMENSION(2)               :: mole_frac_bounds_alpha, mole_frac_bounds_beta, Coa_AB, Caq_AB, Coa_AB_temp, Caq_AB_temp
        REAL(wrp), DIMENSION(2,n_species)     :: activity_coefficient_AB, mass_fraction_water_AB, Coa_j_AB, Caq_j_AB, &
                                                 & Coa_j_AB_temp, Caq_j_AB_temp
        REAL(wrp), DIMENSION(1)               :: mole_a, mole_b, O2C, H2C, Mratio, N2C, func1, func2, ycal_water, ycalc_org_a, &
                                                 & ycalc_org_b, activity_water, activity_calc2_alpha, activity_calc2_beta, mass_fraction_water_a, &
                                                 & mass_fraction_water_b, mass_fraction2, Gibbs_RT, dGibbs_RTdx2, mole_frac_fit, error_out, &
                                                 & guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, aw, C_OA_PM, fit_exit_flag, C_OA_temp, q_alpha_water_temp, &
                                                 & q_alpha_water, guess_C_OA_ugPm3, O2C_temp_2, Mratio_temp_2, H2C_temp_2, aw_point_mean, &
                                                 & weight_q_a, mole_frac_org_a, mole_frac_org_b
        REAL(wrp), DIMENSION(1,1)             :: weight_1_1, aw_1_1, aw_point_1_1
        CHARACTER(clen), DIMENSION(1)         :: BAT_group, BAT_refinement_mode_temp, sim_name, BAT_refinement_mode
        CHARACTER(clen), DIMENSION(n_species) :: BAT_functional_group
        INTEGER(wip)                          :: min_aw_i, s_i, i, j
        REAL(wrp)                             :: avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C
    
        ! Extracting properties
        O2C_values = O2C_j ! O/C ratio of every organic species in the system
        H2C_values = H2C_j ! H/C ratio of every organic species in the system
        N2C_values_densityOnly = N2C_j ! N/C ratio of every organic species in the system
        Molecular_weight = M_j ! Molecular weight of every organic species in the system
        Csat_j_value = Csat_j ! Saturation mass concentration of each organic species in the system
        C_OM_ugPm3 =  Ctotal_j ! Total mass concentration of each organic species in the system
        
        ! Set a value (== 0.0_wrp) for N:C ratio if it is not known at input
        WHERE (N2C_values_densityOnly < -tinynumber) N2C_values_densityOnly = 0.0_wrp
        
        ! Oxygen-bearing functional group of every organic species in the system
        DO i = 1, n_species
            IF (Group_j(i) == 1) THEN
                BAT_functional_group(i) = 'hydroxyl'
            ELSE IF (Group_j(i) == 2) THEN
                BAT_functional_group(i) = 'carboxyl'
            ELSE IF (Group_j(i) == 3) THEN
                BAT_functional_group(i) = 'ketone'
            ELSE IF (Group_j(i) == 4) THEN
                BAT_functional_group(i) = 'hydroperoxide'
            ELSE IF (Group_j(i) == 5) THEN
                BAT_functional_group(i) = 'ether'
            ELSE IF (Group_j(i) == 6) THEN
                BAT_functional_group(i) = 'ester'
            ELSE IF (Group_j(i) == 7) THEN
                BAT_functional_group(i) = 'hydroperoxideSOA'
            ELSE IF (Group_j(i) == 8) THEN
                BAT_functional_group(i) = 'SOA chemicals'
            ELSE IF (Group_j(i) == 9) THEN
                BAT_functional_group(i) = 'PEG'       
            ELSE
                STOP
            END IF
		END DO
		
		! Set simulation options (see SUBROUTINE VBSBAT_setoptions)
        CALL VBSBAT_setoptions(simulation)
        BAT_refinement_mode = simulation%BAT_refinement_mode ! BAT refinement mode
        VBSBAT_options = simulation%VBSBAT_options ! BAT and BAT+VBS simulation options/settings
    
        Molar_mass_ratios = M_water/Molecular_weight ! Molecular weight ratios
    
        ! Main step: check for phase separation conditions and calculation of q_alpha values (fractional liquid-liquid partitioning
        !         of components to phase alpha)
    
        ! Calculate new O/C and Molar mass values using shift fit data to convert functionalized molecules to a hypothetical
        ! OH-equivalent molecule
        CALL convert_chemical_structure_to_OH_eqv_v3(O2C_values, Molar_mass_ratios, BAT_functional_group, O2C_temp, &
            & Mratio_temp)
     
        ! Calculate O/C single phase cross point (O/C at the miscibility limit line as a function of molar mass)
        CALL single_phase_O2C_point_KGv3(Mratio_temp, O2C_single_phase_cross_point_a)
    
        ! Mask: when the converted O/C (OH equiv.) is less than the O/C cross point for the same
        ! molar mass, we are in the less polar region => there is a separation water activity
        O2C_single_phase_cross_point = O2C_values*0.0_wrp
        WHERE (O2C_temp < O2C_single_phase_cross_point_a) O2C_single_phase_cross_point = 1.0_wrp
        mean_prop_mask = O2C_single_phase_cross_point
		
		DO j = 1, n_rh
            aw_series(1) = aw_bat(j) ! Water activity to run the model at
		    aw(1) = aw_series(1) ! Select one water activity walue
            aw_vec = aw(1) ! Create a water activity vector that matches number of organic species
			s_i = 1
    
            IF (VBSBAT_options%force_phase%onePhase == 'no') THEN ! Includes two phase option for q_alpha calculation
                IF (SUM(O2C_single_phase_cross_point) > 0.0_wrp) THEN ! Has a separation water activity separation point
                                                                        ! aw_sep_point because of lower polarity
                    ! Calculate the water activity separation point (a_w_sep_point) for each species
                    ! The a_w_sep_point is determined using the BAT model activities and associated Gibbs energy of mixing
                    ! To approximate the location and aw width over which the liquid–liquid phase separation is prescribed
                    ! to occur, we ﬁrst determine a designated reference point, the so-called water activity separation point
                    ! (aw_sep_point)
                    CALL biphasic_to_single_phase_RH_master_v4(O2C_values, H2C_values, Molar_mass_ratios, &
                        & BAT_functional_group, a_w_sep_point)
                
                    ! Isolate miscible from miscibility gap.
                    a_w_sep_point = a_w_sep_point * O2C_single_phase_cross_point
                
                    ! Creates a a_w_sep_point matrix to match aw_series length and chemical species width.
                    a_w_sep_matrix = SPREAD(a_w_sep_point, 1, 1)
                    aw_series_matrix = SPREAD(aw_series, 2, n_species)
    
                    ! Calculate q_alpha as a function of water activity aw (RH) (sigmoid function to simulate transition between liquid phases)
                    ! The components' fraction in phase alpha follows a smooth transition function for q_alpha with changing RH
                    ! q_alpha_vsRH_values contains q_alpha values of all the species in the system as a function of RH values
                    CALL q_alpha_transfer_vs_aw_calc_v1(a_w_sep_matrix, aw_series_matrix, VBSBAT_options, q_alpha_vsRH_values)
    
                    ! Set threshold: start and end q_alpha values of the sigmoid function
                    max_q_alpha = q_alpha_vsRH_values * 0.0_wrp
                    min_q_alpha = q_alpha_vsRH_values * 0.0_wrp
                    WHERE (q_alpha_vsRH_values > VBSBAT_options%q_alpha%q_alpha_bounds(1)) max_q_alpha = 1.0_wrp
                    WHERE (q_alpha_vsRH_values < VBSBAT_options%q_alpha%q_alpha_bounds(2)) min_q_alpha = 1.0_wrp
    
                    ! Apply q_alpha limits
                    not_min_q_alpha = min_q_alpha * 0.0_wrp
                    WHERE (min_q_alpha < tinynumber) not_min_q_alpha = 1.0_wrp ! where min_q_alpha == 0.0_wrp, not_min_q_alpha = 1.0_wrp
                    q_alpha_vsRH_values = max_q_alpha * q_alpha_vsRH_values
                    q_alpha_vsRH_values = min_q_alpha * q_alpha_vsRH_values + not_min_q_alpha
                ELSE
                    q_alpha_vsRH_values = 1.0 ! Fully miscible, there is no separation water activity
                END IF
    
            ELSE IF (VBSBAT_options%force_phase%onePhase == 'beta') THEN ! Organic-rich phase only (force one phase: q_beta = 1.0_wrp)
                q_alpha_vsRH_values = 0.0_wrp
            ELSE IF (VBSBAT_options%force_phase%onePhase == 'alpha') THEN ! Water-rich phase only (force one phase: q_alpha = 1.0_wrp
                q_alpha_vsRH_values = 1.0_wrp
            ELSE
                PRINT *, 'Specify VBSBAT_options%force_phase%onePhase in SUBROUTINE VBSBAT_options'
            END IF

            ! matrix and vectors to save output data
            partition_coefficients_out = 0.0_wrp
            Cstar_j_out = partition_coefficients_out
            Coa_j_alpha = partition_coefficients_out
            Coa_j_beta = partition_coefficients_out
            Caq_j_alpha = partition_coefficients_out
            Caq_j_beta = partition_coefficients_out
            gamma_alpha = partition_coefficients_out
            gamma_beta = partition_coefficients_out
            mass_fraction_water_alpha_out = partition_coefficients_out
            mass_fraction_water_beta_out = partition_coefficients_out
    
            C_OA_out = 0.0_wrp
            Coa_alpha = C_OA_out
            Coa_beta = C_OA_out
            q_alpha_water_out = C_OA_out
            weight_q_alpha = C_OA_out
            a_w_sep_point_of_meanPM = C_OA_out
    
            ycalc_org_beta_temp = C_OA_out
            mass_fraction_water_beta_temp = C_OA_out
            ycalc_org_alpha_temp = C_OA_out
            mass_fraction_water_alpha_temp = C_OA_out
    
            ! Main step: Use neural network (NN) to get mole fraction of organics at a given water activity (inversion: input: aw ; output: x_org)
            ! We use NN with the BAT model to ﬁnd the correct x_org input for subsequent calculations, since in most applications
            ! aw is known but not x_org
            CALL inverted_NNBAT_v8(O2C_values, H2C_values, Molar_mass_ratios, aw_vec, BAT_functional_group, &
                mole_frac_org_alpha, mole_frac_org_beta)
		    
            mass_fraction_water_alpha = 0.0_wrp
            ycalc_org_alpha = mass_fraction_water_alpha
            mass_fraction_water_beta = 0.0_wrp
            ycalc_org_beta = mass_fraction_water_beta
    
            DO i = 1, n_species ! Start loop through organic species.
            
                ! Set mole fraction limits
                mole_frac_bounds_alpha = [0.0_wrp, 1.0_wrp]
                mole_frac_bounds_beta = [0.0_wrp, 1.0_wrp]
            
                ! IF (VBSBAT_options%BAT_refinement_aw >= aw(1)) THEN
                    ! BAT_refinement_mode_temp = ['none']  ! No interpolation is necessary below RH of BAT_refinement_aw
                ! ELSE
                    ! BAT_refinement_mode_temp = BAT_refinement_mode  ! Interpolate via 501 org mole fraction points when RH is above BAT_refinement_aw
                                                                    ! ! Iterative reﬁnement is required for good agreement with the targeted aw
			    ! END IF
			
			    BAT_refinement_mode_temp = BAT_refinement_mode
                ! Select one species of the system at a time.
                mole_frac_org_a(1) = mole_frac_org_alpha(i) ! Every x_org in phase alpha (water rich)
                mole_frac_org_b(1) = mole_frac_org_beta(i)  ! Every x_org in phase beta (organic-rich)
                O2C(1) = O2C_values(i) ! every O:C
                H2C(1) = H2C_values(i) ! every H:C
                Mratio(1) = Molar_mass_ratios(i) ! every molar mass ratio
                BAT_group(1) = BAT_functional_group(i) ! every funtional group
                N2C(1) = N2C_values_densityOnly(i) ! every N:C
                ycalc_org_a(1) = ycalc_org_alpha(i) ! every org activity coefficient in phase alpha
                ycalc_org_b(1) = ycalc_org_beta(i) ! every org activity coefficient in phase betfractiona
                mass_fraction_water_a(1) = mass_fraction_water_alpha(i) ! every water mass  in phase alpha
                mass_fraction_water_b(1) = mass_fraction_water_beta(i) ! every water mass fraction in phase beta
			
                ! Use x_org to find water mass fraction and activity coefficient associated with each organic species.
                ! Iterative reﬁnement is required for good agreement with the targeted aw (VBSBAT_options%BAT_refinement_aw)
                IF (BAT_refinement_mode_temp(1) == 'interpolate') THEN ! Interpolate in both liquid phases (alpha and beta)
                    BAT_refinement_mode_temp = ['interpolatealpha']
                    CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_a, O2C, H2C, Mratio, BAT_group, &
                        & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_a, &
                        & activity_water, activity_calc2_alpha, mass_fraction_water_a, mass_fraction2, Gibbs_RT, &
                        & dGibbs_RTdx2, mole_frac_fit, error_out)
					
                    BAT_refinement_mode_temp = ['interpolatebeta']
                    CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_b, O2C, H2C, Mratio, BAT_group, &
                        & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_b, &
                        & activity_water, activity_calc2_beta, mass_fraction_water_b, mass_fraction2, Gibbs_RT, &
                        & dGibbs_RTdx2, mole_frac_fit, error_out)
                
			    ELSE ! No refinement is necessary

                    CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_a, O2C, H2C, Mratio, BAT_group, &
                        & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_a, &
                        & activity_water, activity_calc2_alpha, mass_fraction_water_a, mass_fraction2, Gibbs_RT, &
                        & dGibbs_RTdx2, mole_frac_fit, error_out)
                    CALL BAT_activity_calc_with_refinement_v1(mole_frac_org_b, O2C, H2C, Mratio, BAT_group, &
                        & BAT_refinement_mode_temp, aw, N2C, func1, func2, ycal_water, ycalc_org_b, &
                        & activity_water, activity_calc2_beta, mass_fraction_water_b, mass_fraction2, Gibbs_RT, &
                        & dGibbs_RTdx2, mole_frac_fit, error_out)
           
                END IF
                ! Store activity cofficients and water mass fractions (in phase alpha and beta) associated with every
                ! species in the system at a particular RH
                ycalc_org_alpha(i) = ycalc_org_a(1)
                ycalc_org_beta(i) = ycalc_org_b(1)
                mass_fraction_water_alpha(i) = mass_fraction_water_a(1)
                mass_fraction_water_beta(i) = mass_fraction_water_b(1)
            
            END DO ! End loop through organic species
    
            ! y_org_alpha, y_org_beta, w_alpha, w_beta at a given aw
            activity_coefficient_AB(1,:) = ycalc_org_alpha  ! Activity coefficient (of each org) in phase alpha
            activity_coefficient_AB(2,:) = ycalc_org_beta   ! Activity coefficient (of each org) in phase beta
            mass_fraction_water_AB(1,:) = mass_fraction_water_alpha ! Water mass fraction (per org) in phase alpha
            mass_fraction_water_AB(2,:) = mass_fraction_water_beta ! Water mass fraction (per org) in phase beta
    
            ! Main step : VBS
            ! Initial guess for VBS (start point)
            ! Not used if the NN estimates initial guess (set this in VBSBAT_options subroutine)
            guess_C_OAalpha_ugPm3 = 0.0_wrp
            guess_C_OAbeta_ugPm3 = 0.0_wrp
			
			guess_partition_coefficients = 0.0_wrp
			VBSBAT_options%optimization%independent_aw = 'yes'
			VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess = 'yes'
			! Use NN or previous iteration as initial guess
			! ** Check flag first **

			IF (n_rh > 1) THEN
				IF (j > 1) THEN
			        IF (ABS(aw_bat(j)-aw_bat(j-1)) < 0.05_wrp) THEN   ! ** make distinction for different RH domain **
						IF (flag_memory_effect == 1) THEN ! ** Use Logical variable here instead **
				            guess_partition_coefficients = partition_coefficients_memory_effect
						    VBSBAT_options%optimization%independent_aw = 'no'
						END IF
					END IF
				ENDIF
			END IF
	
            ! Single phase calculation (all species in phase beta, sum of q_alpha_org = 0.0_wrp)
            IF (SUM(q_alpha_vsRH_values(s_i,:), 1) < tinynumber) THEN
                q_alpha_values = q_alpha_vsRH_values(s_i,:) ! q_alpha_org of all species in the system at a given RH
    
                ! Solves for partition coefficients Ej (nonideal VBS)
                CALL VBS_equilibration_withLLEpartition_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, &
                    guess_partition_coefficients, C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_values, &
                    mass_fraction_water_AB, Molecular_weight, aw, O2C_values, BAT_functional_group, VBSBAT_options, &
                    partition_coefficients_temp, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, C_OA_PM, q_alpha_water, &
                    error_out)
    
                partition_coefficients = partition_coefficients_temp
                weight_q_alpha(s_i) = q_alpha_vsRH_values(s_i, 1)
    
            ! Single phase calculation (all species in phase alpha, sum of q_alpha_org = number of org species)
            ELSEIF (ABS(SUM(q_alpha_vsRH_values(s_i,:),1)-REAL(n_species, KIND = wrp)) < tinynumber) THEN
                q_alpha_values = q_alpha_vsRH_values(s_i,:)
    
                ! Solves for partition coefficients Ej (nonideal VBS)
                CALL VBS_equilibration_withLLEpartition_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, &
                    guess_partition_coefficients, C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_values, &
                    mass_fraction_water_AB, Molecular_weight, aw, O2C_values, BAT_functional_group, VBSBAT_options, &
                    partition_coefficients_temp, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, C_OA_PM, q_alpha_water, &
                    error_out)
    
                partition_coefficients = partition_coefficients_temp
                weight_q_alpha(s_i) = q_alpha_vsRH_values(s_i, 1)
        
            ELSE ! Calculation for both liquid phases
                ! Beta phase calculation is performed first, which is used for the average q_alpha calculation and
                ! used as a check to see if the two-phase simulation is realistic

                guess_C_OA_ugPm3 = MAX(guess_C_OAbeta_ugPm3, guess_C_OAalpha_ugPm3) ! not used if the NN is used for initial guess

                CALL VBS_equilibration_withLLEpartition_KGv2(guess_C_OA_ugPm3, guess_C_OA_ugPm3, guess_partition_coefficients, &
                    C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_vsRH_values(s_i,:) * 0.0_wrp, mass_fraction_water_AB, &
                    Molecular_weight, aw, O2C_values, BAT_functional_group, VBSBAT_options, partition_coefficients_temp, &
                    Coa_j_AB_temp, Caq_j_AB_temp, Cstar_j_temp, Coa_AB_temp, Caq_AB_temp, C_OA_temp, q_alpha_water_temp, &
                    error_out)
				
                ! Molar based average            
                CALL molar_based_means(Coa_j_AB_temp(2,:), Molecular_weight, O2C_values, H2C_values, N2C_values_densityOnly, &
                & avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C)

                ! Calculate a single water activity separation point based on the average molecular properties for all organics
                CALL biphasic_to_single_phase_RH_master_v4([molar_avg_O2C], [molar_avg_H2C], [M_water/avg_M_gPmol], VBSBAT_options%mean_BAT_functional_group, aw_point_mean)
			
                ! Calculate q_alpha as a function of water activity aw (RH) (sigmoid function to simulate transition between liquid phases)
                aw_point_1_1(1,1) = aw_point_mean(1)
                aw_1_1(1,1) = aw(1)
                CALL q_alpha_transfer_vs_aw_calc_v1(aw_point_1_1, aw_1_1, VBSBAT_options, weight_1_1)
                weight_q_alpha(s_i) = weight_1_1(1,1)
    
                ! Set threshold: start and end q_alpha value for the mean in the sigmoid function
                IF (weight_q_alpha(s_i) > VBSBAT_options%q_alpha%q_alpha_bounds_mean(2)) THEN
                    weight_q_alpha(s_i) = 1.0_wrp
			    END IF
			
                ! Make sure q_alpha mean is greater than mean min, to proceed with 2 phase
                ! calculation as this must pass before a 2 phase calculation is used with
                ! q alpha for each individual species
                IF (weight_q_alpha(s_i) > VBSBAT_options%q_alpha%q_alpha_bounds_mean(1)) THEN
    
                    ! Select which q alpha to use in two phase transition calculation
                    IF (VBSBAT_options%q_alphaVBS%method_to_use == 'mean_prop') THEN
                        q_alpha_molefrac_phase_split_org =  1.0_wrp
                        q_alpha_molefrac_phase_split_org = q_alpha_molefrac_phase_split_org * weight_q_alpha(s_i)
                    ELSE IF (VBSBAT_options%q_alphaVBS%method_to_use == 'individual') THEN
                        q_alpha_molefrac_phase_split_org = q_alpha_vsRH_values(s_i,:)
                    ELSE
                        PRINT *, 'select q_alpha method for transition'
                    END IF
    
                    ! Two phase calculation
                    activity_coefficient_AB(1,:) = ycalc_org_alpha
                    activity_coefficient_AB(2,:) = ycalc_org_beta
                    mass_fraction_water_AB(1,:) = mass_fraction_water_alpha
                    mass_fraction_water_AB(2,:) = mass_fraction_water_beta
    
                    CALL VBS_equilibration_withLLEpartition_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, &
                        guess_partition_coefficients, C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, &
                        q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, Molecular_weight, aw, O2C_values, &
                        BAT_functional_group, VBSBAT_options, partition_coefficients, Coa_j_AB, Caq_j_AB, Cstar_j, &
                        Coa_AB, Caq_AB, C_OA_PM, q_alpha_water, error_out)
            
                    ! Checks if the two phase mass calculation is realistic. Unrealistic cases are identified by the VBS+BAT-predicted liquid
                    ! organic aerosol mass dropping below that predicted for a corresponding single-phase simulation (only beta phase present).
                    IF (VBSBAT_options%q_alphaVBS%Option_to_overwrite_two_phase_calculation == 'yes' &
                        & .AND. (C_OA_PM(1)-C_OA_temp(1))/C_OA_temp(1) < &
                        & VBSBAT_options%q_alphaVBS%overwrite_threshold_two_phase_fraction_difference) THEN
                        ! 2 phase calculation difference is too low
                        ! We use the average of partitioning coefficients Ej of the single-phase prediction and that from a two-phase simulation
                        partition_coefficients = 0.5_wrp * (partition_coefficients + partition_coefficients_temp)
                        Coa_j_AB = 0.5_wrp * (Coa_j_AB_temp + Coa_j_AB)
                        Caq_j_AB = 0.5_wrp * (Caq_j_AB_temp + Caq_j_AB)
                        Cstar_j = 0.5_wrp * (Cstar_j_temp + Cstar_j)
                        Coa_AB = 0.5_wrp * (Coa_AB + Coa_AB_temp)
                        C_OA_PM = 0.5_wrp * (C_OA_PM + C_OA_temp)
                        q_alpha_water = 0.5_wrp * (q_alpha_water + q_alpha_water_temp)
                    END IF

			    ELSE

                    ! Two phase calculation not needed
                    partition_coefficients = partition_coefficients_temp
                    Coa_j_AB = Coa_j_AB_temp
                    Caq_j_AB = Caq_j_AB_temp
                    Cstar_j = Cstar_j_temp
                    Coa_AB = Coa_AB_temp
                    C_OA_PM = C_OA_temp
                    q_alpha_water = q_alpha_water_temp

                END IF

            END IF
    
            ! Store outputs for this aw iteration
            partition_coefficients_out(s_i,:) = partition_coefficients
            Cstar_j_out(s_i,:) = Cstar_j
            Coa_j_alpha(s_i,:) = Coa_j_AB(1,:)
            Coa_j_beta(s_i,:) = Coa_j_AB(2,:)
            Caq_j_alpha(s_i,:) = Caq_j_AB(1,:)
            Caq_j_beta(s_i,:) = Caq_j_AB(2,:)
            gamma_alpha(s_i,:) = ycalc_org_alpha
            gamma_beta(s_i,:) = ycalc_org_beta
            mass_fraction_water_alpha_out(s_i,:) = mass_fraction_water_alpha
            mass_fraction_water_beta_out(s_i,:) = mass_fraction_water_beta
    
            q_alpha_water_out(s_i) = q_alpha_water(1)
            C_OA_out(s_i) = C_OA_PM(1)
            Coa_alpha(s_i) = Coa_AB(1)
            Coa_beta(s_i) = Coa_AB(2)
    
            error_save(s_i) = error_out(1)
            fit_exit_flag_save(s_i) = fit_exit_flag(1)
		
            ! Main step: calculate bulk properties
            Coa_j_PM = Coa_j_alpha + Coa_j_beta ! Total mass concentration of each organic species  in liquid phases (alpha + beta)
            Caq_j_PM = Caq_j_alpha + Caq_j_beta ! Total mass cocentration of water for each organic species in liquid phases (alpha + beta)
    
            C_OA_PM_mat = C_OA_out ! Total mass concentration of organic species in liquid phases (all organics and all liquid phases)
            Caq_PM_mat = SUM(Caq_j_PM, 2) ! Total mass concentration of water in liquid phases (all liquid phases)
    
            ! Outputs
            Coa_org_alpha(j,:) = Coa_j_alpha(1,:) ! Mass concentration of each organics in phases alpha and beta (ug/m3)
            Coa_org_beta(j,:) =  Coa_j_beta(1,:)
            Coa_total_water_alpha(j) = SUM(Caq_j_alpha) ! Total mass concentration of water in phases alpha and beta (ug/m3)
            Coa_total_water_beta(j) = SUM(Caq_j_beta)
			
		END DO ! end water activity loop

	END SUBROUTINE BAT_INVERSION_VBS_SIMULATION
		
    !****************************************************************************************
    !*   SUBROUTINE VBS_equilibration_withLLEpartition_KGv2     							*
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Solves for VBS partition coefficients Ej                                           *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************

    SUBROUTINE VBS_equilibration_withLLEpartition_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3,  &
        & guess_partition_coefficients, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, &
        & mass_fraction_water_AB, molecular_weight, a_water, O2C_values, BAT_functional_group, VBSBAT_options, &
        & partition_coefficients_AB, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: guess_C_OAalpha_ugPm3, & ! Guess total organic mass concentration in phase alpha
                                                     & guess_C_OAbeta_ugPm3, & ! Guess total organic mass concentration in phase beta
                                                     & guess_partition_coefficients, & ! Guess partitioning coefficient of each org. species
                                                     & C_OM_ugPm3, & ! Total mass concentration of organics (particle+gas)
                                                     & Cstar_dry, & ! Saturation mass concentration of each org. species
                                                     & q_alpha_molefrac_phase_split_org, & ! Fractional liquid-liquid partitioning of org. species to                                   ! phase alpha
                                                     & molecular_weight, & ! Molecular weight of each org. species
                                                     & O2C_values, & ! Elemental O/C ratio of each org. species
                                                     & a_water ! Water activity
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: BAT_functional_group ! Main oxygen-bearing functional group of each organic species
        TYPE(VBSBAT_options_input), INTENT(IN)    :: VBSBAT_options ! Simulation options
        REAL(wrp), DIMENSION(:,:), INTENT(IN)     :: activity_coefficient_AB, & ! Activity coefficients of organic species in phases alpha and beta
                                                     & mass_fraction_water_AB ! Water mass fraction associated with each org. species in liquid
                                                                              ! phases alpha and beta

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: partition_coefficients_AB, & ! Partitioning coefficients of organics Ej
                                                     & Cstar_j, & ! Effective saturation concentration of organics C*j
                                                     & Coa, & ! Total organic mass concentration in the particle phases
                                                     & q_alpha_water, & ! Fractional liquid-liquid partitioning of water to phase alpha
                                                     & Coa_AB, & ! Total organic mass concentration in phases alpha and beta
                                                     & Caq_AB, & ! Total water mass concentration in phases alpha and beta
                                                     & error_out !Eerror for convergence of VBS
        REAL(wrp), DIMENSION(:,:), INTENT(OUT)    :: Coa_j_AB, & ! Organic mass concentration of each organic species in the particle phase
                                                     & Caq_j_AB ! Water mass concentration associated with each organic species in the particle phase

        ! Local variables
        REAL(wrp) :: guess_C_oaaqbeta_ugPm3
        REAL(wrp), DIMENSION(SIZE(molecular_weight,1)) :: Molar_mass_ratios, Cstar_j_guess, O2C_temp, Mratio_temp, Ej_guess, &
                                                          & fit_CoaEj
        REAL(wrp), DIMENSION(1)                        :: guess_error

        ! KIND of the variables below is predefined. Used in hybrd1 (third party code)
        REAL(8), DIMENSION(SIZE(molecular_weight,1))   :: Ej_guess8, fvec
        REAL(8)                                        :: newerror_out8 , guess_error8, tol
        INTEGER(4)                                     :: info, n, prin, i
    
        ! Molecular weight ratios (Mwater/Morg)
        Molar_mass_ratios = M_water/molecular_weight
    
        ! Get a guess for the total mass concentration in phase beta (water + organic): guess_C_oaaqbeta_ugPm3
        IF (guess_C_OAbeta_ugPm3(1) > 0.0_wrp + tinynumber) THEN
            guess_C_oaaqbeta_ugPm3 = guess_C_OAbeta_ugPm3(1)/(1.0_wrp-SUM(mass_fraction_water_AB(2,:), 1) / &
                & SIZE(molecular_weight,1))
        ELSE IF (guess_C_OAalpha_ugPm3(1) > 0.0_wrp + tinynumber) THEN
            guess_C_oaaqbeta_ugPm3 = guess_C_OAalpha_ugPm3(1)/(1.0_wrp-SUM(mass_fraction_water_AB(2,:),1) / &
                & SIZE(molecular_weight,1))
        ELSE
            guess_C_oaaqbeta_ugPm3 = 1.0_wrp/(1.0_wrp-SUM(mass_fraction_water_AB(2,:),1)/SIZE(molecular_weight,1))
		END IF
		
  !      ! Estimate Ej guess: guess partitioning coefficients for organic species j.
  !      IF (SUM(guess_partition_coefficients, 1) < 0.0_wrp + tinynumber) THEN
  !          IF (VBSBAT_options%optimization%independent_aw == 'yes') THEN ! Do not use guess for Ej from previous water activity iteration
  !
  !              IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'no') THEN ! Estimate guess for Ej using an approximation
  !                  Cstar_j_guess = Cstar_dry * activity_coefficient_AB(2,:)/(1.0_wrp + mass_fraction_water_AB(2,:)*(1.0_wrp / &
  !                      & Molar_mass_ratios - 1.0_wrp))
  !                  Ej_guess = 1.0_wrp/(1.0_wrp + Cstar_j_guess/guess_C_oaaqbeta_ugPm3)
  !
		!		ELSE IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'yes')  THEN ! Guess for Ej using Neural Network
  !                  ! Network trained on hydroxyl so need to convert to OH first
  !                  CALL convert_chemical_structure_to_OH_eqv_v3(O2C_values, Molar_mass_ratios, BAT_functional_group, O2C_temp, & !bpt here
  !                      Mratio_temp)
  !                  CALL VBSBAT_neural_network_estimate_v1(Cstar_dry, C_OM_ugPm3, O2C_temp, M_water/Mratio_temp, &
  !                      mass_fraction_water_AB(2,:), a_water, VBSBAT_options%VBSBAT_NN_options, Ej_guess) ! correct neural network
  !              END IF
  !          END IF
		!ELSE IF (VBSBAT_options%optimization%independent_aw == 'no') THEN ! Use previous Ej_guess results
  !          Ej_guess = guess_partition_coefficients 
  !      ELSE
  !          PRINT *, 'Select VBSBAT_options%optimization%independent_aw option: either yes or no'
		!END IF
		
        ! Estimate Ej guess: guess partitioning coefficients for organic species j.

		IF (VBSBAT_options%optimization%independent_aw == 'yes') THEN ! Do not use guess for Ej from previous water activity iteration
            IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'no') THEN ! Estimate guess for Ej using an approximation
                Cstar_j_guess = Cstar_dry * activity_coefficient_AB(2,:)/(1.0_wrp + mass_fraction_water_AB(2,:)*(1.0_wrp / &
                    & Molar_mass_ratios - 1.0_wrp))
                Ej_guess = 1.0_wrp/(1.0_wrp + Cstar_j_guess/guess_C_oaaqbeta_ugPm3)

			ELSE IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'yes')  THEN ! Guess for Ej using Neural Network
                ! Network trained on hydroxyl so need to convert to OH first
                CALL convert_chemical_structure_to_OH_eqv_v3(O2C_values, Molar_mass_ratios, BAT_functional_group, O2C_temp, & !bpt here
                    Mratio_temp)
                CALL VBSBAT_neural_network_estimate_v1(Cstar_dry, C_OM_ugPm3, O2C_temp, M_water/Mratio_temp, &
                    mass_fraction_water_AB(2,:), a_water, VBSBAT_options%VBSBAT_NN_options, Ej_guess) ! correct neural network

            END IF

		ELSE IF (VBSBAT_options%optimization%independent_aw == 'no') THEN ! Use previous Ej_guess results
		    Ej_guess = guess_partition_coefficients 

		END IF
		
        ALLOCATE(C_OM_ugPm3_objfun(SIZE(molecular_weight,1)), Cstar_dry_objfun(SIZE(molecular_weight,1)), activity_coefficient_AB_objfun(2,SIZE(molecular_weight,1)), &
            & mass_fraction_water_AB_objfun(2,SIZE(molecular_weight,1)), q_alpha_molefrac_phase_split_org_objfun(SIZE(molecular_weight,1)), molecular_weight_objfun(SIZE(molecular_weight,1)))
        ALLOCATE(partition_coefficients_objfun(SIZE(molecular_weight,1)), Cstar_j_objfun(SIZE(molecular_weight,1)), Coa_AB_objfun(SIZE(molecular_weight,1)), &
			& Caq_AB_objfun(SIZE(molecular_weight,1)), Coa_objfun(SIZE(molecular_weight,1)), q_alpha_water_objfun(SIZE(molecular_weight,1)), Coa_j_AB_objfun(2,SIZE(molecular_weight,1)), &
			& Caq_j_AB_objfun(2,SIZE(molecular_weight,1)))

        ! Assign variables declared inside MODULE OBJECTIVE_FUNCTION_MOD
        C_OM_ugPm3_objfun = C_OM_ugPm3 ! Total mass concentration of organics (particle+gas)
        Cstar_dry_objfun = Cstar_dry ! Saturation mass concentration of each chemical species
        activity_coefficient_AB_objfun = activity_coefficient_AB ! Activity coefficient of each chemical species in liquid phases alpha and beta
        mass_fraction_water_AB_objfun = mass_fraction_water_AB ! Mass fraction of the water associated with each organic species in phases alpha and beta
        q_alpha_molefrac_phase_split_org_objfun = q_alpha_molefrac_phase_split_org ! Fractional liquid-liquid partitioning of org. to phase alpha
        molecular_weight_objfun = molecular_weight ! Molecular weight of each chemical species

        ! Assign 0.0_wrp to variables
        partition_coefficients_objfun = 0.0_wrp * C_OM_ugPm3 ! Partitioning coefficients Ej
        Cstar_j_objfun = 0.0_wrp * C_OM_ugPm3 ! Effective saturation mass concentration of each chemical species C*j
        Coa_AB_objfun = 0.0_wrp * C_OM_ugPm3 ! Total mass concentration of organics in liquid phases alpha and beta
        Caq_AB_objfun = 0.0_wrp * C_OM_ugPm3 ! Total mass concentration of water in liquid phases alpha and beta
        Coa_objfun  = 0.0_wrp ! Total organic mass concentration (phase alpha + phase beta)
        q_alpha_water_objfun = 0.0_wrp ! Fractional liquid-liquid partitioning of water to phase alpha for each org. j
        Coa_j_AB_objfun = 0.0_wrp * activity_coefficient_AB ! Mass concentration of organics species j in liquid phases alpha and beta
        Caq_j_AB_objfun = 0.0_wrp * activity_coefficient_AB ! Mass concentration of water (for each org. j) in liquid phases alpha and beta
        
        ! Calculate error in guess
        CALL VBS_equilibration_withLLEpartition_objFun_KGv3(Ej_guess, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
            & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight, partition_coefficients_AB, &
            & Coa_j_AB, Caq_j_AB, Cstar_j,  Coa_AB, Caq_AB, Coa, q_alpha_water, guess_error)
        
        ! IF (VBSBAT_options%optimization%guess_refinement_threshold <= guess_error(1)) THEN  ! Check if refinement is needed
            ! !IF (VBSBAT_options%optimization%opt_method == 'powell') THEN ! Use Powell's hybrid method (hybrd1) since refinement is needed
            ! !    n = SIZE(molecular_weight)
            ! !    Ej_guess8 = REAL(Ej_guess, KIND=8)
            ! !    tol = REAL(VBSBAT_options%optimization%fit_tolerance, KIND = 8)
            ! !    CALL hybrd1(error_out_powell, n, Ej_guess8, fvec, tol, info)
            ! !    fit_CoaEj = REAL(Ej_guess8, KIND = wrp)
            ! !ELSE IF (VBSBAT_options%optimization%opt_method == 'none') THEN !  Use Ej_guess directly
            ! !    fit_CoaEj = Ej_guess
            ! !ELSE
            ! !    WRITE(*,*) 'Invalid VBSBAT_options%optimization%opt_method. options in VBSBAT_options or settings input file: none or powell'
            ! !END IF

		IF (VBSBAT_options%optimization%opt_method == 'powell') THEN ! Use Powell's hybrid method (hybrd1) since refinement is needed
			n = SIZE(molecular_weight)
			tol = VBSBAT_options%optimization%fit_tolerance
			CALL hybrd1(error_out_powell, n, Ej_guess, fvec, tol, info)
			fit_CoaEj = Ej_guess
			flag_memory_effect = info   ! ** have a variable that counts the frequency of failures 
			partition_coefficients_memory_effect = Ej_guess
		ELSE IF (VBSBAT_options%optimization%opt_method == 'none') THEN !  Use Ej_guess directly
			fit_CoaEj = Ej_guess
		ELSE
			WRITE(*,*) 'Invalid VBSBAT_options%optimization%opt_method. options in VBSBAT_options or settings input file: none or powell'
		END IF

        ! ELSE
            ! fit_CoaEj = Ej_guess ! No refinement is neeeded
        ! END IF

        ! Runs a final calculation for output
        CALL VBS_equilibration_withLLEpartition_objFun_KGv3(fit_CoaEj, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
            & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight, partition_coefficients_AB, &
            Coa_j_AB, Caq_j_AB, Cstar_j,  Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        DEALLOCATE(C_OM_ugPm3_objfun, Cstar_dry_objfun, activity_coefficient_AB_objfun, &
            mass_fraction_water_AB_objfun, q_alpha_molefrac_phase_split_org_objfun, molecular_weight_objfun)
        DEALLOCATE(partition_coefficients_objfun, Cstar_j_objfun, Coa_AB_objfun, Caq_AB_objfun, Coa_objfun, &
            q_alpha_water_objfun, Coa_j_AB_objfun, Caq_j_AB_objfun)

    END SUBROUTINE VBS_equilibration_withLLEpartition_KGv2

!****************************************************************************************
!*   SUBROUTINE VBS_equilibration_withLLEpartition_objFun_KGv3					        *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Calculates the error for convergence of VBS, based on squared error of Ej.         *
!*   Calculates a new Ej passed on a guess Ej, and outputs the phase masses.            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 06-05-2021                                                      *
!****************************************************************************************

    SUBROUTINE VBS_equilibration_withLLEpartition_objFun_KGv3(guess_Ej, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
        & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight,  partition_coefficients, Coa_j_AB, &
        & Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)    :: guess_Ej, & ! Guess partitioning coefficient Ej of each org. species
                                                  & C_OM_ugPm3, & ! Total mass concentration of organics (particle+gas)
                                                  & Cstar_dry, & ! Saturation mass concentration of each org. species
                                                  & q_alpha_molefrac_phase_split_org, & ! Fractional liquid-liquid partitioning of org. species to
                                                                                        ! phase alpha
                                                  & molecular_weight ! Molecular weight of each org. species
        REAL(wrp), DIMENSION(:,:), INTENT(IN)  :: mass_fraction_water_AB, & ! Water mass fraction associated with each org. species in liquid phases
                                                                            ! alpha and beta
                                                  & activity_coefficient_AB ! Activity coefficients of organic species in phases alpha and beta

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: partition_coefficients, & ! Partitioning coefficients of organics Ej in phases alpha and beta
                                                     & Cstar_j, & ! Effective saturation concentration of organics C*j
                                                     & Coa, & ! Total organic mass concentration in the particle phases (alpha+beta)
                                                     & q_alpha_water, & ! Fractional liquid-liquid partitioning of water to phase alpha
                                                     & Coa_AB, & ! Total organic mass concentration in phases alpha and beta
                                                     & Caq_AB, & ! Total water mass concentration in phases alpha and beta
                                                     & error_out !Eerror for convergence of VBS
        REAL(wrp), DIMENSION(:,:), INTENT(OUT)    :: Coa_j_AB, & ! Organic mass concentration of each organic species in the particle phase
                                                     & Caq_j_AB ! Water mass concentration associated with each organic species in the particle phase

        ! Local variables
        REAL(wrp), DIMENSION(SIZE(molecular_weight, 1)) :: mass_fraction_water_alpha, mass_fraction_water_beta, activity_coefficient_alpha, &
                                                           & activity_coefficient_beta, q_alpha, q_beta, mass_fraction_water_alpha_denominator, &
                                                           & mass_fraction_water_beta_denominator, alpha_good_denominators, &
                                                           & not_alpha_good_denominators, beta_good_denominators, not_beta_good_denominators, &
                                                           & Coa_j, Coa_j_alpha, Coa_j_beta, Cstar_j_via_alpha, Cstar_j_via_beta, Cstar_j_used, &
                                                           & Ej_new, Coa_j_new, Coa_j_alpha_new, Caq_j_alpha_new, Coa_j_beta_new, Caq_j_beta_new
        REAL(wrp), DIMENSION(1)                         :: Coa_guess_viaEj, Caq_alpha, Coaaq_alpha, Caq_beta, Coaaq_beta, Coaaq_via_Ej_guess, &
                                                           & massweighted_molar_weight_alpha, massweighted_molar_weight_beta, Coa_new_viaEj, Caq_alpha_new, &
                                                           & Coa_alpha_new, Caq_beta_new, Coa_beta_new, q_alpha_water_new

        ! Get terms for alpha and beta
        mass_fraction_water_alpha = mass_fraction_water_AB(1,:) ! Water mass fraction alpha and beta (all species) for a given aw
        mass_fraction_water_beta = mass_fraction_water_AB(2,:)
        activity_coefficient_alpha = activity_coefficient_AB(1,:) ! Activiy coefficients alpha and beta (all species) for a give aw
        activity_coefficient_beta = activity_coefficient_AB(2,:)

        q_alpha = q_alpha_molefrac_phase_split_org
        q_beta = 1.0_wrp - q_alpha_molefrac_phase_split_org

        ! checks that denomintor is not Inf or NaN
        mass_fraction_water_alpha_denominator = (1.0_wrp - mass_fraction_water_alpha)
        mass_fraction_water_beta_denominator = (1.0_wrp - mass_fraction_water_beta)
        
        WHERE (mass_fraction_water_alpha_denominator < tinynumber) mass_fraction_water_alpha_denominator = 1.0_wrp
        WHERE (mass_fraction_water_beta_denominator < tinynumber) mass_fraction_water_beta_denominator = 1.0_wrp

        ! New Coa j calculation
        Coa_j = guess_Ej * C_OM_ugPm3
        Coa_guess_viaEj(1) = SUM(Coa_j, 1)

        ! alpha phase calculation
        Coa_j_alpha = Coa_j * q_alpha ! organic mass
        Caq_alpha(1) = SUM(Coa_j_alpha * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator, 1) ! water mass
        ! based on organic mass and mass fraction of water from BAT mode
        Coaaq_alpha(1) = SUM(Coa_j_alpha, 1) + Caq_alpha(1) ! water plus organic mass
        
        ! beta phase calculation
        Coa_j_beta = Coa_j * q_beta ! organic mass
        Caq_beta(1) = SUM(Coa_j_beta * mass_fraction_water_beta / mass_fraction_water_beta_denominator, 1) ! water mass
        ! based on organic mass and mass fraction of water from BAT mode
        Coaaq_beta(1) = SUM(Coa_j_beta, 1) + Caq_beta(1) ! water plus organic mass
       
        Coaaq_via_Ej_guess(1) = Coaaq_beta(1) + Coaaq_alpha(1) ! total organic plus water mass

		
        ! C* via alpha phase
        massweighted_molar_weight_alpha(1) = SUM((Coa_j_alpha/(molecular_weight)), 1) + Caq_alpha(1)/M_water

        IF (massweighted_molar_weight_alpha(1) > 0.0_wrp + tinynumber) THEN
            Cstar_j_via_alpha = Cstar_dry * activity_coefficient_alpha * q_alpha * Coaaq_via_Ej_guess(1) / &
			    & (molecular_weight * massweighted_molar_weight_alpha(1))
        ELSE
            Cstar_j_via_alpha = q_alpha * 0.0_wrp
        END IF

        ! C* via beta
        massweighted_molar_weight_beta(1) = SUM((Coa_j_beta/(molecular_weight)), 1) + Caq_beta(1)/M_water

        IF (massweighted_molar_weight_beta(1) > 0.0_wrp) THEN
            Cstar_j_via_beta = Cstar_dry * activity_coefficient_beta * q_beta * Coaaq_via_Ej_guess(1) / &
			    & (molecular_weight * massweighted_molar_weight_beta(1))
        ELSE
            Cstar_j_via_beta = q_beta * 0.0_wrp
        END IF

        ! select C* to use
        ! could have different methods used here
        ! weight average by q_alpha
        Cstar_j_used = Cstar_j_via_alpha * q_alpha + Cstar_j_via_beta * q_beta

        ! new Ej calculation
        Ej_new = (1.0_wrp + Cstar_j_used/Coaaq_via_Ej_guess(1))**(-1)

        ! calculate new mass values
        Coa_j_new = Ej_new * C_OM_ugPm3
        Coa_new_viaEj(1) = SUM(Coa_j_new, 1)

        ! alpha masses
        Coa_j_alpha_new = Coa_j * q_alpha
        Caq_j_alpha_new = (Coa_j_alpha_new * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator)
        Caq_alpha_new(1) = SUM(Caq_j_alpha_new, 1)
        Coa_alpha_new(1) = SUM(Coa_j_alpha_new, 1)

        ! beta masses
        Coa_j_beta_new = Coa_j_new * q_beta
        Caq_j_beta_new = (Coa_j_beta_new * mass_fraction_water_beta / mass_fraction_water_beta_denominator)
        Caq_beta_new(1) = SUM(Caq_j_beta_new, 1)
        Coa_beta_new(1) = SUM(Coa_j_beta_new, 1)

        ! Calculate the new fraction of water for each organic species j in liquid phase alpha
        IF (Caq_alpha_new(1) + Caq_beta_new(1) < tinynumber) THEN
            q_alpha_water_new(1) = 0.0_wrp
        ELSE
            q_alpha_water_new(1) = Caq_alpha_new(1)/(Caq_alpha_new(1) + Caq_beta_new(1))
		END IF
	
        ! error for convergence of VBS, based on squared error of Ej and organic ! mass
        error_out(1) = SUM(((guess_Ej - Ej_new)**2 + (Coa_guess_viaEj(1) - Coa_new_viaEj(1))**2),1)

        ! collect the outputs
        partition_coefficients = Ej_new
        q_alpha_water = q_alpha_water_new

        Coa_j_AB(1,:) = Coa_j_alpha_new
        Coa_j_AB(2,:) = Coa_j_beta_new
        Caq_j_AB(1,:) = Caq_j_alpha_new
        Caq_j_AB(2,:) = Caq_j_beta_new

        Coa_AB(1) = Coa_alpha_new(1)
        Coa_AB(2) = Coa_beta_new(1)
        Caq_AB(1) = Caq_alpha_new(1)
        Caq_AB(2) = Caq_beta_new(1)

        Coa(1) = SUM(Coa_AB, 1)
        Cstar_j = Cstar_j_used

	END SUBROUTINE VBS_equilibration_withLLEpartition_objFun_KGv3
		
    !****************************************************************************************
    !*   SUBROUTINE CSAT_ESTIMATION					                                        *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Estimates the pure component saturation mass concentration of an organic species   *
	!*	 j, C_sat_j (ug/m3) at the temperature of interest, at dry conditions (RH = 0),     *
	!*   when C_star_j (ug/m3) and C_PM_j (ug/m3) are known.								*
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 19-02-2022                                                             *
    !*   -> latest changes: 23-02-2022                                                      *
    !****************************************************************************************
    
    SUBROUTINE CSAT_ESTIMATION(n_species, M_j, O2C_j, Cstar_j, Ctotal_j, Group_j, Csat_j)
    
        IMPLICIT NONE
    
        ! Input
        INTEGER(wip), INTENT(IN)                        :: n_species                        ! Number of organic species in the system
        INTEGER(wip), DIMENSION(n_species), INTENT(IN)  :: Group_j                          ! Type of oxygen-bearing functional group of each org. species
        REAL(wrp), DIMENSION(n_species), INTENT(IN)     :: M_j, O2C_j, &					! Molar mass, O/C elemental ratios of each org. species
                                                          & Cstar_j, Ctotal_j				! Total mass concentration (gas+particle) and saturation mass
                                                                                            ! concentration of each org. species       
        
        ! Output      
		REAL(wrp), DIMENSION(n_species), INTENT(OUT)    :: Csat_j                           ! Pure componenent saturation concentration at the temperature of interest (ug/m3)
		
		! Local variables
        REAL(wrp)       							    :: Coa_total_water_alpha, Coa_total_water_beta
        REAL(wrp), DIMENSION(n_species)				    :: Coa_org_alpha, Coa_org_beta 
        TYPE(input_data)							    :: simulation
        TYPE(specialoptions_input)                      :: special_options
        TYPE(growth_input)                              :: growth
        TYPE(VBSBAT_options_input)                      :: VBSBAT_options
        REAL(wrp), DIMENSION(1)                         :: C_OA_PM_mat, Caq_PM_mat
        REAL(wrp), DIMENSION(n_species)                 :: N2C_values_densityOnly, Molar_mass_ratios, O2C_eqv, molarmass_ratio_eqv, &
                                                            & O2C_single_phase_cross_point_a, O2C_single_phase_cross_point, mean_prop_mask, a_w_sep_point, aw_vec, &
                                                            & mole_frac_org_alpha, mole_frac_org_beta, mass_fraction_water_alpha, ycalc_org_alpha, ycalc_org_beta, &
                                                            & mass_fraction_water_beta, guess_partition_coefficients, partition_coefficients_temp, Cstar_j_temp, &
                                                            & partition_coefficients, mass_inPM_temp, mass_fraction_inPM, q_alpha_molefrac_phase_split_org, &
                                                            & Mratio_temp, O2C_temp, q_alpha_values, Cstar_j_values, C_OM_ugPm3, O2C_values, H2C_values, Molar_mass
        REAL(wrp), DIMENSION(1)                         :: aw_series, C_OA_out, Coa_alpha, Coa_beta, q_alpha_water_out, weight_q_alpha, &
                                                            & a_w_sep_point_of_meanPM, ycalc_org_beta_temp, mass_fraction_water_beta_temp, ycalc_org_alpha_temp, &
                                                            & mass_fraction_water_alpha_temp, fit_exit_flag_save, error_save, C_OA_ratio, kappa
        REAL(wrp), DIMENSION(1,n_species)               :: a_w_sep_matrix, aw_series_matrix, q_alpha_vsRH_values, max_q_alpha, &
                                                            & min_q_alpha, not_min_q_alpha, partition_coefficients_out, Cstar_j_out, Coa_j_alpha, Coa_j_beta, Caq_j_alpha, &
                                                            & Caq_j_beta, gamma_alpha, gamma_beta, mass_fraction_water_alpha_out, mass_fraction_water_beta_out, Coa_j_PM, &
                                                            & Caq_j_PM
        REAL(wrp), DIMENSION(n_species+1)               :: Molecular_weight_withwater
        REAL(wrp), DIMENSION(2)                         :: mole_frac_bounds_alpha, mole_frac_bounds_beta, Coa_AB, Caq_AB, Coa_AB_temp, Caq_AB_temp
        REAL(wrp), DIMENSION(2,n_species)               :: activity_coefficient_AB, mass_fraction_water_AB, Coa_j_AB, Caq_j_AB, &
                                                            & Coa_j_AB_temp, Caq_j_AB_temp
        REAL(wrp), DIMENSION(1)                         :: mole_a, mole_b, O2C, H2C, Mratio, N2C, func1, func2, ycal_water, ycalc_org_a, &
                                                            & ycalc_org_b, activity_water, activity_calc2_alpha, activity_calc2_beta, mass_fraction_water_a, &
                                                            & mass_fraction_water_b, mass_fraction2, Gibbs_RT, dGibbs_RTdx2, mole_frac_fit, error_out, &
                                                            & guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, aw, C_OA_PM, fit_exit_flag, C_OA_temp, q_alpha_water_temp, &
                                                            & q_alpha_water, guess_C_OA_ugPm3, O2C_temp_2, Mratio_temp_2, H2C_temp_2, aw_point_mean, &
                                                            & weight_q_a, mole_frac_org_a, mole_frac_org_b
        REAL(wrp), DIMENSION(1,1)                       :: weight_1_1, aw_1_1, aw_point_1_1
        CHARACTER(clen), DIMENSION(1)                   :: BAT_group, BAT_refinement_mode_temp, sim_name, BAT_refinement_mode
        CHARACTER(clen), DIMENSION(n_species)           :: BAT_functional_group
        INTEGER(wip)                                    :: min_aw_i, s_i, i
        REAL(wrp)                                       :: avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C
    
        ! Extracting properties
        ! Oxygen-bearing functional group of every organic species in the system
        DO i = 1, n_species
            IF (Group_j(i) == 1) THEN
                BAT_functional_group(i) = 'hydroxyl'
            ELSE IF (Group_j(i) == 2) THEN
                BAT_functional_group(i) = 'carboxyl'
            ELSE IF (Group_j(i) == 3) THEN
                BAT_functional_group(i) = 'ketone'
            ELSE IF (Group_j(i) == 4) THEN
                BAT_functional_group(i) = 'hydroperoxide'
            ELSE IF (Group_j(i) == 5) THEN
                BAT_functional_group(i) = 'ether'
            ELSE IF (Group_j(i) == 6) THEN
                BAT_functional_group(i) = 'ester'
            ELSE IF (Group_j(i) == 7) THEN
                BAT_functional_group(i) = 'hydroperoxideSOA'
            ELSE IF (Group_j(i) == 8) THEN
                BAT_functional_group(i) = 'SOA chemicals'
            ELSE IF (Group_j(i) == 9) THEN
                BAT_functional_group(i) = 'PEG'       
            ELSE
                STOP
            END IF
        END DO
        
        ! Set simulation options (see SUBROUTINE VBSBAT_setoptions)
        CALL VBSBAT_setoptions(simulation)
        BAT_refinement_mode = simulation%BAT_refinement_mode ! BAT refinement mode
        VBSBAT_options = simulation%VBSBAT_options ! BAT and BAT+VBS simulation options/settings
		
		aw(1) = 0.0_wrp ! Water activity to run the model at
		
        ! y_org_alpha, y_org_beta, w_alpha, w_beta at a given aw
        activity_coefficient_AB(1,:) = 0.0_wrp  ! Activity coefficient (of each org) in phase alpha
        activity_coefficient_AB(2,:) = 1.0_wrp   ! Activity coefficient (of each org) in phase beta
        mass_fraction_water_AB(1,:) = 0.0_wrp ! Water mass fraction (per org) in phase alpha
        mass_fraction_water_AB(2,:) = 0.0_wrp ! Water mass fraction (per org) in phase beta
    
        ! Main step : VBS
        ! Initial guess for VBS (start point)
        ! Not used if the NN estimates initial guess (set this in VBSBAT_options subroutine)
        guess_C_OAalpha_ugPm3 = 0.1_wrp
        guess_C_OAbeta_ugPm3 = 0.1_wrp
        guess_partition_coefficients = 0.0_wrp
        
        ! Single phase calculation (all species in phase beta, sum of q_alpha_org = 0.0_wrp)
		 q_alpha_values = 0.0_wrp ! q_alpha_org of all species in the system at a given RH
    
		! Solves for partition coefficients Ej (nonideal VBS)
        CALL VBS_equilibration_withLLEpartition_Csat_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3, &
            guess_partition_coefficients, Ctotal_j, Cstar_j, activity_coefficient_AB, q_alpha_values, &
            mass_fraction_water_AB, M_j, aw, O2C_j, BAT_functional_group, VBSBAT_options, &
            partition_coefficients_temp, Coa_j_AB, Caq_j_AB, Cstar_j_temp, Coa_AB, Caq_AB, C_OA_PM, q_alpha_water, &
            error_out)
    
        partition_coefficients = partition_coefficients_temp
		
		! Calculate an approximate value for Csat at dry conditions
		CALL VBS_equilibration_extractCsat_withLLEpartition_KGv2(Coa_j_AB(2,:), Cstar_j, &
		    & M_j,  Csat_j)
		
    END SUBROUTINE CSAT_ESTIMATION

    !****************************************************************************************
    !*   SUBROUTINE VBS_equilibration_withLLEpartition_Csat_KGv2							*
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Solves for VBS partition coefficients Ej                                           *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************

    SUBROUTINE VBS_equilibration_withLLEpartition_Csat_KGv2(guess_C_OAalpha_ugPm3, guess_C_OAbeta_ugPm3,  &
        & guess_partition_coefficients, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, &
        & mass_fraction_water_AB, molecular_weight, a_water, O2C_values, BAT_functional_group, VBSBAT_options, &
        & partition_coefficients_AB, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: guess_C_OAalpha_ugPm3, & ! Guess total organic mass concentration in phase alpha
                                                     & guess_C_OAbeta_ugPm3, & ! Guess total organic mass concentration in phase beta
                                                     & guess_partition_coefficients, & ! Guess partitioning coefficient of each org. species
                                                     & C_OM_ugPm3, & ! Total mass concentration of organics (particle+gas)
                                                     & Cstar_dry, & ! Saturation mass concentration of each org. species
                                                     & q_alpha_molefrac_phase_split_org, & ! Fractional liquid-liquid partitioning of org. species to                                   ! phase alpha
                                                     & molecular_weight, & ! Molecular weight of each org. species
                                                     & O2C_values, & ! Elemental O/C ratio of each org. species
                                                     & a_water ! Water activity
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: BAT_functional_group ! Main oxygen-bearing functional group of each organic species
        TYPE(VBSBAT_options_input), INTENT(IN)    :: VBSBAT_options ! Simulation options
        REAL(wrp), DIMENSION(:,:), INTENT(IN)     :: activity_coefficient_AB, & ! Activity coefficients of organic species in phases alpha and beta
                                                     & mass_fraction_water_AB ! Water mass fraction associated with each org. species in liquid
                                                                              ! phases alpha and beta

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: partition_coefficients_AB, & ! Partitioning coefficients of organics Ej
                                                     & Cstar_j, & ! Effective saturation concentration of organics C*j
                                                     & Coa, & ! Total organic mass concentration in the particle phases
                                                     & q_alpha_water, & ! Fractional liquid-liquid partitioning of water to phase alpha
                                                     & Coa_AB, & ! Total organic mass concentration in phases alpha and beta
                                                     & Caq_AB, & ! Total water mass concentration in phases alpha and beta
                                                     & error_out !Eerror for convergence of VBS
        REAL(wrp), DIMENSION(:,:), INTENT(OUT)    :: Coa_j_AB, & ! Organic mass concentration of each organic species in the particle phase
                                                     & Caq_j_AB ! Water mass concentration associated with each organic species in the particle phase

        ! Local variables
        REAL(wrp) :: guess_C_oaaqbeta_ugPm3
        REAL(wrp), DIMENSION(SIZE(molecular_weight,1)) :: Molar_mass_ratios, Cstar_j_guess, O2C_temp, Mratio_temp, Ej_guess, &
                                                          & fit_CoaEj
        REAL(wrp), DIMENSION(1)                        :: guess_error

        ! KIND of the variables below is predefined. Used in hybrd1 (third party code)
        REAL(8), DIMENSION(SIZE(molecular_weight,1))   :: Ej_guess8, fvec
        REAL(8)                                        :: newerror_out8 , guess_error8, tol
        INTEGER(4)                                     :: info, n, prin, i
		    
        ! Molecular weight ratios (Mwater/Morg)
        Molar_mass_ratios = M_water/molecular_weight
    
        ! Get a guess for the total mass concentration in phase beta (water + organic): guess_C_oaaqbeta_ugPm3
        IF (guess_C_OAbeta_ugPm3(1) > 0.0_wrp + tinynumber) THEN
            guess_C_oaaqbeta_ugPm3 = guess_C_OAbeta_ugPm3(1)/(1.0_wrp-SUM(mass_fraction_water_AB(2,:), 1) / &
                & SIZE(molecular_weight,1))
        ELSE IF (guess_C_OAalpha_ugPm3(1) > 0.0_wrp + tinynumber) THEN
            guess_C_oaaqbeta_ugPm3 = guess_C_OAalpha_ugPm3(1)/(1.0_wrp-SUM(mass_fraction_water_AB(2,:),1) / &
                & SIZE(molecular_weight,1))
        ELSE
            guess_C_oaaqbeta_ugPm3 = 1.0_wrp/(1.0_wrp-SUM(mass_fraction_water_AB(2,:),1)/SIZE(molecular_weight,1))
		END IF
		
		

        ! Estimate Ej guess: guess partitioning coefficients for organic species j.
        IF (SUM(guess_partition_coefficients, 1) < 0.0_wrp + tinynumber) THEN
            IF (VBSBAT_options%optimization%independent_aw == 'yes') THEN ! Do not use guess for Ej from previous water activity iteration

                IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'no') THEN ! Estimate guess for Ej using an approximation
                    Cstar_j_guess = Cstar_dry * activity_coefficient_AB(2,:)/(1.0_wrp + mass_fraction_water_AB(2,:)*(1.0_wrp / &
                        & Molar_mass_ratios - 1.0_wrp))
                    Ej_guess = 1.0_wrp/(1.0_wrp + Cstar_j_guess/guess_C_oaaqbeta_ugPm3)

				ELSE IF (VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess == 'yes')  THEN ! Guess for Ej using Neural Network
                    ! Network trained on hydroxyl so need to convert to OH first
                    CALL convert_chemical_structure_to_OH_eqv_v3(O2C_values, Molar_mass_ratios, BAT_functional_group, O2C_temp, & !bpt here
                        Mratio_temp)
                    CALL VBSBAT_neural_network_estimate_v1(Cstar_dry, C_OM_ugPm3, O2C_temp, M_water/Mratio_temp, &
                        mass_fraction_water_AB(2,:), a_water, VBSBAT_options%VBSBAT_NN_options, Ej_guess) ! correct neural network
                END IF
            END IF
		ELSE IF (VBSBAT_options%optimization%independent_aw == 'no') THEN ! Use previous Ej_guess results
            Ej_guess = guess_partition_coefficients 
        ELSE
            PRINT *, 'Select VBSBAT_options%optimization%independent_aw option: either yes or no'
		END IF
		
        ALLOCATE(C_OM_ugPm3_objfun(SIZE(molecular_weight,1)), Cstar_dry_objfun(SIZE(molecular_weight,1)), activity_coefficient_AB_objfun(2,SIZE(molecular_weight,1)), &
            & mass_fraction_water_AB_objfun(2,SIZE(molecular_weight,1)), q_alpha_molefrac_phase_split_org_objfun(SIZE(molecular_weight,1)), molecular_weight_objfun(SIZE(molecular_weight,1)))
        ALLOCATE(partition_coefficients_objfun(SIZE(molecular_weight,1)), Cstar_j_objfun(SIZE(molecular_weight,1)), Coa_AB_objfun(SIZE(molecular_weight,1)), &
			& Caq_AB_objfun(SIZE(molecular_weight,1)), Coa_objfun(SIZE(molecular_weight,1)), q_alpha_water_objfun(SIZE(molecular_weight,1)), Coa_j_AB_objfun(2,SIZE(molecular_weight,1)), &
			& Caq_j_AB_objfun(2,SIZE(molecular_weight,1)))

        ! Assign variables declared inside MODULE OBJECTIVE_FUNCTION_MOD
        C_OM_ugPm3_objfun = C_OM_ugPm3 ! Total mass concentration of organics (particle+gas)
        Cstar_dry_objfun = Cstar_dry ! Saturation mass concentration of each chemical species
        activity_coefficient_AB_objfun = activity_coefficient_AB ! Activity coefficient of each chemical species in liquid phases alpha and beta
        mass_fraction_water_AB_objfun = mass_fraction_water_AB ! Mass fraction of the water associated with each organic species in phases alpha and beta
        q_alpha_molefrac_phase_split_org_objfun = q_alpha_molefrac_phase_split_org ! Fractional liquid-liquid partitioning of org. to phase alpha
        molecular_weight_objfun = molecular_weight ! Molecular weight of each chemical species

        ! Assign 0.0_wrp to variables
        partition_coefficients_objfun = 0.0_wrp * C_OM_ugPm3 ! Partitioning coefficients Ej
        Cstar_j_objfun = 0.0_wrp * C_OM_ugPm3 ! Effective saturation mass concentration of each chemical species C*j
        Coa_AB_objfun = 0.0_wrp * C_OM_ugPm3 ! Total mass concentration of organics in liquid phases alpha and beta
        Caq_AB_objfun = 0.0_wrp * C_OM_ugPm3 ! Total mass concentration of water in liquid phases alpha and beta
        Coa_objfun  = 0.0_wrp ! Total organic mass concentration (phase alpha + phase beta)
        q_alpha_water_objfun = 0.0_wrp ! Fractional liquid-liquid partitioning of water to phase alpha for each org. j
        Coa_j_AB_objfun = 0.0_wrp * activity_coefficient_AB ! Mass concentration of organics species j in liquid phases alpha and beta
        Caq_j_AB_objfun = 0.0_wrp * activity_coefficient_AB ! Mass concentration of water (for each org. j) in liquid phases alpha and beta
        
        ! ! Calculate error in guess
        ! CALL VBS_equilibration_withLLEpartition_objFun_Csat_KGv3(Ej_guess, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
            ! & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight, partition_coefficients_AB, &
            ! & Coa_j_AB, Caq_j_AB, Cstar_j,  Coa_AB, Caq_AB, Coa, q_alpha_water, guess_error)
        
        ! IF (VBSBAT_options%optimization%guess_refinement_threshold <= guess_error(1)) THEN  ! Check if refinement is needed
            ! !IF (VBSBAT_options%optimization%opt_method == 'powell') THEN ! Use Powell's hybrid method (hybrd1) since refinement is needed
            ! !    n = SIZE(molecular_weight)
            ! !    Ej_guess8 = REAL(Ej_guess, KIND=8)
            ! !    tol = REAL(VBSBAT_options%optimization%fit_tolerance, KIND = 8)
            ! !    CALL hybrd1(error_out_powell, n, Ej_guess8, fvec, tol, info)
            ! !    fit_CoaEj = REAL(Ej_guess8, KIND = wrp)
            ! !ELSE IF (VBSBAT_options%optimization%opt_method == 'none') THEN !  Use Ej_guess directly
            ! !    fit_CoaEj = Ej_guess
            ! !ELSE
            ! !    WRITE(*,*) 'Invalid VBSBAT_options%optimization%opt_method. options in VBSBAT_options or settings input file: none or powell'
            ! !END IF

		
		IF (VBSBAT_options%optimization%opt_method == 'powell') THEN ! Use Powell's hybrid method (hybrd1) since refinement is needed
			n = SIZE(molecular_weight)
			tol = VBSBAT_options%optimization%fit_tolerance
			CALL hybrd1(error_out_powell_Csat, n, Ej_guess, fvec, tol, info)
			!PRINT *, info
			fit_CoaEj = Ej_guess
		ELSE IF (VBSBAT_options%optimization%opt_method == 'none') THEN !  Use Ej_guess directly
			fit_CoaEj = Ej_guess
		ELSE
			WRITE(*,*) 'Invalid VBSBAT_options%optimization%opt_method. options in VBSBAT_options or settings input file: none or powell'
		END IF
        
        ! ELSE
            ! fit_CoaEj = Ej_guess ! No refinement is neeeded
        ! END IF

        ! Runs a final calculation for output
        CALL VBS_equilibration_withLLEpartition_objFun_Csat_KGv3(fit_CoaEj, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
            & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight, partition_coefficients_AB, &
            Coa_j_AB, Caq_j_AB, Cstar_j,  Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        DEALLOCATE(C_OM_ugPm3_objfun, Cstar_dry_objfun, activity_coefficient_AB_objfun, &
            mass_fraction_water_AB_objfun, q_alpha_molefrac_phase_split_org_objfun, molecular_weight_objfun)
        DEALLOCATE(partition_coefficients_objfun, Cstar_j_objfun, Coa_AB_objfun, Caq_AB_objfun, Coa_objfun, &
            q_alpha_water_objfun, Coa_j_AB_objfun, Caq_j_AB_objfun)

    END SUBROUTINE VBS_equilibration_withLLEpartition_Csat_KGv2

!****************************************************************************************
!*   SUBROUTINE VBS_equilibration_withLLEpartition_objFun_Csat_KGv3					    *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Calculates the error for convergence of VBS, based on squared error of Ej.         *
!*   Calculates a new Ej passed on a guess Ej, and outputs the phase masses.            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 06-05-2021                                                      *
!****************************************************************************************

    SUBROUTINE VBS_equilibration_withLLEpartition_objFun_Csat_KGv3(guess_Ej, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, &
        & q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight,  partition_coefficients, Coa_j_AB, &
        & Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)    :: guess_Ej, & ! Guess partitioning coefficient Ej of each org. species
                                                  & C_OM_ugPm3, & ! Total mass concentration of organics (particle+gas)
                                                  & Cstar_dry, & ! Saturation mass concentration of each org. species
                                                  & q_alpha_molefrac_phase_split_org, & ! Fractional liquid-liquid partitioning of org. species to
                                                                                        ! phase alpha
                                                  & molecular_weight ! Molecular weight of each org. species
        REAL(wrp), DIMENSION(:,:), INTENT(IN)  :: mass_fraction_water_AB, & ! Water mass fraction associated with each org. species in liquid phases
                                                                            ! alpha and beta
                                                  & activity_coefficient_AB ! Activity coefficients of organic species in phases alpha and beta

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: partition_coefficients, & ! Partitioning coefficients of organics Ej in phases alpha and beta
                                                     & Cstar_j, & ! Effective saturation concentration of organics C*j
                                                     & Coa, & ! Total organic mass concentration in the particle phases (alpha+beta)
                                                     & q_alpha_water, & ! Fractional liquid-liquid partitioning of water to phase alpha
                                                     & Coa_AB, & ! Total organic mass concentration in phases alpha and beta
                                                     & Caq_AB, & ! Total water mass concentration in phases alpha and beta
                                                     & error_out !Eerror for convergence of VBS
        REAL(wrp), DIMENSION(:,:), INTENT(OUT)    :: Coa_j_AB, & ! Organic mass concentration of each organic species in the particle phase
                                                     & Caq_j_AB ! Water mass concentration associated with each organic species in the particle phase

        ! Local variables
        REAL(wrp), DIMENSION(SIZE(molecular_weight, 1)) :: mass_fraction_water_alpha, mass_fraction_water_beta, activity_coefficient_alpha, &
                                                           & activity_coefficient_beta, q_alpha, q_beta, mass_fraction_water_alpha_denominator, &
                                                           & mass_fraction_water_beta_denominator, alpha_good_denominators, &
                                                           & not_alpha_good_denominators, beta_good_denominators, not_beta_good_denominators, &
                                                           & Coa_j, Coa_j_alpha, Coa_j_beta, Cstar_j_via_alpha, Cstar_j_via_beta, Cstar_j_used, &
                                                           & Ej_new, Coa_j_new, Coa_j_alpha_new, Caq_j_alpha_new, Coa_j_beta_new, Caq_j_beta_new
        REAL(wrp), DIMENSION(1)                         :: Coa_guess_viaEj, Caq_alpha, Coaaq_alpha, Caq_beta, Coaaq_beta, Coaaq_via_Ej_guess, &
                                                           & massweighted_molar_weight_alpha, massweighted_molar_weight_beta, Coa_new_viaEj, Caq_alpha_new, &
                                                           & Coa_alpha_new, Caq_beta_new, Coa_beta_new, q_alpha_water_new

        ! ! Get terms for alpha and beta
        ! mass_fraction_water_alpha = mass_fraction_water_AB(1,:) ! Water mass fraction alpha and beta (all species) for a given aw
        ! mass_fraction_water_beta = mass_fraction_water_AB(2,:)
        ! activity_coefficient_alpha = activity_coefficient_AB(1,:) ! Activiy coefficients alpha and beta (all species) for a give aw
        ! activity_coefficient_beta = activity_coefficient_AB(2,:)

        ! q_alpha = q_alpha_molefrac_phase_split_org
        ! q_beta = 1.0_wrp - q_alpha_molefrac_phase_split_org

        ! ! checks that denomintor is not Inf or NaN
        ! mass_fraction_water_alpha_denominator = (1.0_wrp - mass_fraction_water_alpha)
        ! mass_fraction_water_beta_denominator = (1.0_wrp - mass_fraction_water_beta)
        
        ! WHERE (mass_fraction_water_alpha_denominator < tinynumber) mass_fraction_water_alpha_denominator = 1.0_wrp
        ! WHERE (mass_fraction_water_beta_denominator < tinynumber) mass_fraction_water_beta_denominator = 1.0_wrp

        ! New Coa j calculation
        Coa_j = guess_Ej * C_OM_ugPm3
        Coa_guess_viaEj(1) = SUM(Coa_j, 1)

        ! ! alpha phase calculation
        ! Coa_j_alpha = Coa_j * q_alpha ! organic mass
        ! Caq_alpha(1) = SUM(Coa_j_alpha * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator, 1) ! water mass
        ! ! based on organic mass and mass fraction of water from BAT mode
        ! Coaaq_alpha(1) = SUM(Coa_j_alpha, 1) + Caq_alpha(1) ! water plus organic mass
        
        ! ! beta phase calculation
        ! Coa_j_beta = Coa_j * q_beta ! organic mass
        ! Caq_beta(1) = SUM(Coa_j_beta * mass_fraction_water_beta / mass_fraction_water_beta_denominator, 1) ! water mass
        ! ! based on organic mass and mass fraction of water from BAT mode
        ! Coaaq_beta(1) = SUM(Coa_j_beta, 1) + Caq_beta(1) ! water plus organic mass
       
        !Coaaq_via_Ej_guess(1) = Coaaq_beta(1) + Coaaq_alpha(1) ! total organic plus water mass
		Coaaq_via_Ej_guess(1) = Coa_guess_viaEj(1)
		
        ! ! C* via alpha phase
        ! massweighted_molar_weight_alpha(1) = SUM((Coa_j_alpha/(molecular_weight)), 1) + Caq_alpha(1)/M_water

        ! IF (massweighted_molar_weight_alpha(1) > 0.0_wrp + tinynumber) THEN
            ! Cstar_j_via_alpha = Cstar_dry * activity_coefficient_alpha * q_alpha !* Coaaq_via_Ej_guess(1) / &
			    ! !& (molecular_weight * massweighted_molar_weight_alpha(1))
        ! ELSE
            ! Cstar_j_via_alpha = q_alpha * 0.0_wrp
        ! END IF

        ! ! C* via beta
        ! massweighted_molar_weight_beta(1) = SUM((Coa_j_beta/(molecular_weight)), 1) + Caq_beta(1)/M_water

        ! IF (massweighted_molar_weight_beta(1) > 0.0_wrp) THEN
            ! Cstar_j_via_beta = Cstar_dry * activity_coefficient_beta * q_beta !* Coaaq_via_Ej_guess(1) / &
			    ! !& (molecular_weight * massweighted_molar_weight_beta(1))
        ! ELSE
            ! Cstar_j_via_beta = q_beta * 0.0_wrp
        ! END IF
		
        ! select C* to use
        ! could have different methods used here
        ! weight average by q_alpha
        !Cstar_j_used = Cstar_j_via_alpha * q_alpha + Cstar_j_via_beta * q_beta
		
		Cstar_j_used = Cstar_dry

        ! new Ej calculation
        Ej_new = (1.0_wrp + Cstar_j_used/Coaaq_via_Ej_guess(1))**(-1)

        ! calculate new mass values
        Coa_j_new = Ej_new * C_OM_ugPm3
        Coa_new_viaEj(1) = SUM(Coa_j_new, 1)

        ! ! alpha masses
        ! Coa_j_alpha_new = Coa_j * q_alpha
        ! Caq_j_alpha_new = (Coa_j_alpha_new * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator)
        ! Caq_alpha_new(1) = SUM(Caq_j_alpha_new, 1)
        ! Coa_alpha_new(1) = SUM(Coa_j_alpha_new, 1)
	    Coa_j_alpha_new = 0.0_wrp
        Caq_j_alpha_new = 0.0_wrp
        Caq_alpha_new(1) = 0.0_wrp
        Coa_alpha_new(1) = 0.0_wrp

        ! ! beta masses
        ! Coa_j_beta_new = Coa_j_new * q_beta
        ! Caq_j_beta_new = (Coa_j_beta_new * mass_fraction_water_beta / mass_fraction_water_beta_denominator)
        ! Caq_beta_new(1) = SUM(Caq_j_beta_new, 1)
        ! Coa_beta_new(1) = SUM(Coa_j_beta_new, 1)
        Coa_j_beta_new = Coa_j_new
        Caq_j_beta_new = 0.0_wrp
        Caq_beta_new(1) = 0.0_wrp
        Coa_beta_new(1) = Coa_new_viaEj(1)

        ! Calculate the new fraction of water for each organic species j in liquid phase alpha
        IF (Caq_alpha_new(1) + Caq_beta_new(1) < tinynumber) THEN
            q_alpha_water_new(1) = 0.0_wrp
        ELSE
            q_alpha_water_new(1) = Caq_alpha_new(1)/(Caq_alpha_new(1) + Caq_beta_new(1))
		END IF
	
        ! error for convergence of VBS, based on squared error of Ej and organic mass		
        error_out(1) = SUM(((guess_Ej - Ej_new)**2 + (Coa_guess_viaEj(1) - Coa_new_viaEj(1))**2),1)
		
        ! collect the outputs
        partition_coefficients = Ej_new
        q_alpha_water = q_alpha_water_new

        Coa_j_AB(1,:) = Coa_j_alpha_new
        Coa_j_AB(2,:) = Coa_j_beta_new
        Caq_j_AB(1,:) = Caq_j_alpha_new
        Caq_j_AB(2,:) = Caq_j_beta_new

        Coa_AB(1) = Coa_alpha_new(1)
        Coa_AB(2) = Coa_beta_new(1)
        Caq_AB(1) = Caq_alpha_new(1)
        Caq_AB(2) = Caq_beta_new(1)

        Coa(1) = Coa_new_viaEj(1)
        Cstar_j = Cstar_j_used

	END SUBROUTINE VBS_equilibration_withLLEpartition_objFun_Csat_KGv3

    !****************************************************************************************
    !*   VBS_equilibration_extractCsat_withLLEpartition_KGv2							    *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Estimates the pure component saturation mass concentration of an organic species   *
	!*	 j, C_sat_j (ug/m3) at the temperature of interest, at dry conditions (RH = 0),     *
	!*   when C_star_j (ug/m3) and C_PM_j (ug/m3) are known.								*
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 19-02-2022                                                             *
    !*   -> latest changes: 23-02-2022                                                      *
    !****************************************************************************************
      
    SUBROUTINE VBS_equilibration_extractCsat_withLLEpartition_KGv2(C_PM_org_j, C_star_org_j, &
		& M_org_j,  C_sat_org_j_approx)

	IMPLICIT NONE
	
	    ! Inputs
	    REAL(wrp), DIMENSION(:), INTENT(IN)             ::  C_PM_org_j, C_star_org_j, M_org_j

	    ! Outputs
	    REAL(wrp), DIMENSION(:), INTENT(OUT)			::  C_sat_org_j_approx
        ! Locals
		REAL(wrp)										::  C_PM_org_water
		
		! In this subroutine, we are using equation (7) of Gorkowski et al., 2019 to estimate C_sat_j from C_star_j and C_PM_org_j at dry conditions; q_beta_org_j = 1, y_beta_org_j = 1
		
		! At dry conditions, there is not water in the particle and only an organic-rich phase exists: C_PM_org_water = C_PM_beta_org
		C_PM_org_water = SUM(C_PM_org_j)
		
		C_sat_org_j_approx = C_star_org_j/(C_PM_org_water/(M_org_j * SUM(C_PM_org_j/M_org_j)))
	
	END SUBROUTINE VBS_equilibration_extractCsat_withLLEpartition_KGv2
	

    !****************************************************************************************
    !*   SUBROUTINE VBSBAT_setoptions                                                       *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Sets the different options/settings for BAT and BAT+VBS simulations.               *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes: 06-05-2021                                                      *
    !****************************************************************************************
      
    SUBROUTINE VBSBAT_setoptions(simulation)

        IMPLICIT NONE

        ! Input and Output
        TYPE(input_data), INTENT(INOUT) :: simulation ! Variable that stores the options for the BAT and BAT+VBS simulations

        ! Default options
        ALLOCATE(simulation%BAT_refinement_mode(1), simulation%VBSBAT_options%run_mode_used(1))

        ! BAT options
        simulation%VBSBAT_options%BAT_refinement_aw = 0.90_wrp ! Water activity above which an iterative
                                                               ! reﬁnement is required for good agreement with the
                                                               ! targeted water activity.
                                                               ! Needed due to errors in neural network fit.
                                                               ! Iterative reﬁnement of the mole fraction of organics
                                                               ! to match the given water activity.
        simulation%BAT_refinement_mode(1) = 'interpolate' ! Use interpolation to determine the mole fraction of organics above a
                                                          ! relative humidity of BAT_refinement_aw.

        simulation%VBSBAT_options%mean_BAT_functional_group='hydroxyl' ! Average functional group used to calculate the water activity
                                                                       ! separation point of the multi-organic mixture.
                                                                       ! Options: 'hydroxyl', 'carboxyl', 'ether', 'ketone', 'ester',
                                                                       ! 'hydroperoxide', 'hydroperoxideSOA', 'PEG'

        ! q_alpha options
        simulation%VBSBAT_options%force_phase%onePhase = 'no' ! phase calculation isolation (q calculation)
                                                              ! force one phase in q calculaiton or not
                                                              ! options = 'alpha','beta','no'
        simulation%VBSBAT_options%q_alpha%min_spread_in_aw = 1.0E-06_wrp ! Min water activity gap from the separation point aw to complete aqueous
                                                                         ! dilution (where aw -> 1)
        simulation%VBSBAT_options%q_alpha%q_alpha_at_1phase_aw = 0.99_wrp ! Value of q_alpha_org (= q_alpha_water_sep) at the water activity
                                                                          ! separation point
        simulation%VBSBAT_options%q_alpha%q_alpha_bounds = [1.0E-06_wrp, 1.0_wrp] ! Sets q_alpha limits
        simulation%VBSBAT_options%q_alpha%q_alpha_bounds_mean = [1.0E-06_wrp, 1.0_wrp] ! Set q alpha limits for mean PM
        simulation%VBSBAT_options%q_alphaVBS%method_to_use = 'individual' ! select high q alpha to use in two phase
                                                                          ! transition calculation
                                                                          ! options = 'mean_prop', 'individual'

        ! Used to check if the two phase mass calculation is realistic
        simulation%VBSBAT_options%q_alphaVBS%Option_to_overwrite_two_phase_calculation = 'yes' ! Overwrite if calculation is unrelistic
                                                                                               ! options = 'yes', 'no'
        simulation%VBSBAT_options%q_alphaVBS%overwrite_threshold_two_phase_fraction_difference = 0.0_wrp ! threshold to compare 2-phase vs 1-phase
                                                                                                         ! simulations
        ! VBS
        simulation%VBSBAT_options%optimization%independent_aw = 'yes'  ! Use guess partitioning that is dependent on lower water activity iteration
                                                                       ! This is only relevant when using an array of aw as input (MATLAB version)
                                                                       ! options = 'yes' or 'no'
        simulation%VBSBAT_options%VBSBAT_NN_options%use_NN_for_VBS_initial_guess = 'yes' ! options = 'yes' or 'no'
        simulation%VBSBAT_options%VBSBAT_NN_options%NN_type = 'individual_properties' ! options = 'mean_properties', 'individual_properties'
		simulation%VBSBAT_options%optimization%guess_refinement_threshold = SQRT(EPSILON(1.0_wrp)) ! refinement criteria of NN guess for partitioning coefficients
        simulation%VBSBAT_options%optimization%opt_method = 'powell' ! solver options = 'none', 'powell'
        simulation%VBSBAT_options%optimization%fit_tolerance = SQRT(EPSILON(1.0_wrp))  ! solver tolerance

        simulation%VBSBAT_options%run_mode_used(1) = 'default' ! The Python interface uses the default option

        ! Refinement only when RH > 90 %
        IF (simulation%VBSBAT_options%run_mode_used(1) == 'default') THEN
        ! No change

        ! Refinement at any RH; the tolerance of the solver is lower; 
        ELSE IF (simulation%VBSBAT_options%run_mode_used(1) == 'robust') THEN
            simulation%VBSBAT_options%BAT_refinement_aw = 0.0_wrp
			simulation%VBSBAT_options%optimization%guess_refinement_threshold = SQRT(EPSILON(1.0_wrp))
			simulation%VBSBAT_options%optimization%opt_method = 'powell'
            simulation%VBSBAT_options%optimization%fit_tolerance = SQRT(EPSILON(1.0_wrp))
		END IF

    END SUBROUTINE VBSBAT_setoptions

END MODULE BATVBS_MOD
