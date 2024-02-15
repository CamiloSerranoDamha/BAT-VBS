!****************************************************************************************
!*   MODULE DATATYPE_MOD                                                                *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Derived data type                                                                  *
!*   Module that define derived data types used for BAT and BAT+VBS simulations         *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************
    
MODULE DATATYPE_MOD

    USE PRECISION_MOD ! Module that defines the precision of REAL, INTEGER data types and length of CHAR

    IMPLICIT NONE

    TYPE :: q_alpha_input
        REAL(wrp) :: min_spread_in_aw, q_alpha_at_1phase_aw
        REAL(wrp), DIMENSION(2) :: q_alpha_bounds, q_alpha_bounds_mean
    END TYPE q_alpha_input

    TYPE :: q_alphaVBS_input
        CHARACTER(clen) :: Option_to_overwrite_two_phase_calculation, method_to_use
        CHARACTER(clen), DIMENSION(2) :: method_to_use_options
        REAL(wrp) :: overwrite_threshold_two_phase_fraction_difference
    END TYPE q_alphaVBS_input

    TYPE :: force_phase_input
        CHARACTER(clen) :: onePhase
        CHARACTER(clen), DIMENSION(3) :: onePhase_options
    END TYPE force_phase_input

    TYPE :: optimization_input
        REAL(wrp)       :: fit_tolerance, guess_refinement_threshold
        INTEGER              :: MaxIter
        CHARACTER(clen)      :: opt_method, independent_aw
    END TYPE optimization_input

    TYPE :: VBSBAT_NN_options_input
        CHARACTER(clen) :: use_NN_for_VBS_initial_guess, NN_type
    END TYPE VBSBAT_NN_options_input

    TYPE :: VBSBAT_options_input
        TYPE(q_alpha_input) :: q_alpha
        TYPE(q_alphaVBS_input) :: q_alphaVBS
        TYPE(force_phase_input) :: force_phase
        TYPE(optimization_input) :: optimization
        TYPE(VBSBAT_NN_options_input) :: VBSBAT_NN_options
        CHARACTER(clen), DIMENSION(8) :: BAT_functional_group_options
        CHARACTER(clen), DIMENSION(5) :: run_mode_possible_options
        CHARACTER(clen), DIMENSION(:), ALLOCATABLE :: plot_PM, run_mode_used
        CHARACTER(clen), DIMENSION(1) :: mean_BAT_functional_group
        REAL(wrp) :: BAT_refinement_aw, BAT_refinement_tolerance
    END TYPE VBSBAT_options_input

    TYPE :: system_input
    REAL(wrp), DIMENSION(:), ALLOCATABLE :: calculate_Csat_j_with_aw
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: Molecular_weight, O2C_values, H2C_values, Csat_j_value, &
            C_OM_ugPm3, optional_Cstar_ugPm3, optional_Cliquid_ugPm3
        CHARACTER(clen), DIMENSION(:), ALLOCATABLE :: BAT_functional_group, speciesName
    END TYPE system_input

    TYPE :: input_data
        CHARACTER(clen), DIMENSION(:), ALLOCATABLE :: BAT_refinement_mode, run_name
        TYPE(VBSBAT_options_input) :: VBSBAT_options
        TYPE(system_input) :: system
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: water_activity
    END TYPE input_data


    TYPE :: kappa_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: kappaHGF, kappa
        CHARACTER(clen) :: note
    END TYPE kappa_input

    TYPE :: Organic_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: volume, radius, mass, mean_density_ugPm3
    END TYPE Organic_input

    TYPE :: PM_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: volume_HGF, radius_HGF, mass_HGF, volume, radius, mass
    END TYPE PM_input

    TYPE :: settings_input
        REAL(wrp) :: Diameter_vol_eqv_org_nm, sigma_droplet_N_m, density_water_g_cm3, Mw, R, Temp, sigmaW_droplet_N_m, &
            sigmaOrg_droplet_N_m, fixed_aw_kappaCCN
        CHARACTER(clen) :: surface_tension_method
    END TYPE settings_input

    TYPE :: beta_input
        REAL(wrp) :: phase
    END TYPE beta_input

    TYPE :: alpha_input
        REAL(wrp) :: aw, SatCritical, kappa_SatCritical, kappa_a_w, Dp_critical, V_w_critical
    END TYPE alpha_input

    TYPE :: fixed_aw_input
        REAL(wrp) :: fixed_aw_val
        TYPE(beta_input) :: beta
        TYPE(alpha_input) :: alpha
    END TYPE fixed_aw_input

    TYPE :: kappaCCN_direct_input
        TYPE(beta_input) :: beta
        TYPE(alpha_input) :: alpha
        TYPE(fixed_aw_input) :: fixed_aw
        TYPE(settings_input) :: settings
    END TYPE kappaCCN_direct_input

    TYPE :: growth_input
        TYPE(kappa_input) :: kappa
        TYPE(Organic_input) :: Organic
        TYPE(PM_input) :: PM
        TYPE(kappaCCN_direct_input) :: kappaCCN_direct
    END TYPE growth_input

    TYPE :: postEquilb_input
        REAL(wrp), DIMENSION(:,:), ALLOCATABLE :: mole_fraction_water_free_alpha, mole_fraction_water_free_beta, &
            mole_fraction_alpha, mole_fraction_beta
    END TYPE postEquilb_input

    TYPE :: species_specific_input
        REAL(wrp), DIMENSION(:,:), ALLOCATABLE :: partition_coefficients, Cstar_j, Coa_j_alpha, Coa_j_beta, &
            Caq_j_alpha, Caq_j_beta, Coa_j_PM, Caq_j_PM, gamma_alpha, gamma_beta, mass_fraction_water_alpha, &
            mass_fraction_water_beta, q_alpha_org_values
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: q_alpha_water, q_alpha_org_used, mass_fraction_inPM
        TYPE(postEquilb_input) :: postEquilb
    END TYPE species_specific_input

    TYPE :: totals_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: Coa_alpha, Coa_beta, C_OA_PM, Caq_PM, Caq_PM_alpha, Caq_PM_beta, &
            C_OA_ratio, C_OA
    END TYPE totals_input

    TYPE :: kappas_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: kappaHGF, kappa
    END TYPE kappas_input

    TYPE :: inputs_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: Csat_j_value, C_OM_ugPm3, O2C_values, H2C_values, Molecular_weight
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: aw_series
        CHARACTER(clen), DIMENSION(:), ALLOCATABLE :: BAT_functional_group
        CHARACTER(clen) :: BAT_refinement_mode
        TYPE(VBSBAT_options_input) :: VBSBAT_options
    END TYPE inputs_input

    TYPE :: details_input
        TYPE(species_specific_input) :: species_specific
        TYPE(totals_input) :: totals
        TYPE(growth_input) :: growth
        TYPE(kappas_input) :: kappas
        TYPE(inputs_input) :: inputs
    END TYPE details_input

    TYPE :: specialoptions_input
        CHARACTER(clen) :: fit_option
        REAL(wrp) :: tran_lowO2C_fractionOne_phase, tran_lowO2C_sigmoid_bend, tran_lowO2C_sigmoid_shift, &
            tran_midO2C_sigmoid_bend, tran_midO2C_sigmoid_shift
    END TYPE specialoptions_input

    TYPE :: step_input
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: gain
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xoffset
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: ymin
    END TYPE step_input

    TYPE :: kappa_settings_input
        REAL(wrp) :: Diameter_vol_eqv_org_nm, sigma_droplet_N_m, density_water_g_cm3, Mw, R, Temp, sigmaW_droplet_N_m, &
            sigmaOrg_droplet_N_m, fixed_aw_kappaCCN
        CHARACTER(clen) :: surface_tension_method
    END TYPE kappa_settings_input

END MODULE DATATYPE_MOD


