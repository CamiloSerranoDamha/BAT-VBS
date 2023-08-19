!****************************************************************************************
!*   MODULE BAT_MOD                                                                     *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that includes all the subroutines related to the prediction of the extent   *
!*   of liquid-liquid phase separation, as well as the calculation of activity          *
!*  coefficients and mass fractions of water and organics.                              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 20-10-2022                                                      *
!****************************************************************************************

MODULE BAT_MOD

    USE PRECISION_MOD ! Module that defines the precision/length of REAL/INTERGER/CHAR data type
    USE DATATYPE_MOD  ! Module that defines derived data types
    USE TOOLS_MOD     ! Module that contains subroutines to perform different recurrent calculations

    IMPLICIT NONE

CONTAINS

    !****************************************************************************************
    !*   SUBROUTINE biphasic_to_single_phase_RH_master_v4                                   *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Finds the water activity at which liquid-liquid phase separation takes place.      *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE biphasic_to_single_phase_RH_master_v4(O2C, H2C, Mratio, BAT_functional_group, RH_cross_point)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: O2C, H2C, Mratio ! Elemental ratios (O/C, H/C) and molecular weight ratios (Mwater/Morg) of
                                                                      ! organic species
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: BAT_functional_group ! Representative oxygen-bearing functinal group of each organic species

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: RH_cross_point ! Relative humidity over which liquid-liquid phase separation happens

        ! Local variables
        REAL(wrp), DIMENSION(1)                       :: phase_sep_check
        REAL(wrp), DIMENSION(SIZE(RH_cross_point, 1)) :: lower_round, not_lower_round, upper_round, not_upper_round
        INTEGER(wip), DIMENSION(1)                    :: lower_a_w_sep_index, upper_a_w_sep_index, matching_Upper_a_w_sep_index
        REAL(wrp), DIMENSION(501)                     :: func1, func2, ycalc_water, ycalc_org, activity_water, activity_org, &
                                                         & mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2, mole_frac
        INTEGER(wip)                                  :: i, interpolate_step_numb

        RH_cross_point = 0.0_wrp

        ! Number of interpolation points
        interpolate_step_numb = 500_wip

        ! Create 501 mole fraction points separated by 1/500
        mole_frac(1) = 0.0_wrp
        mole_frac(SIZE(mole_frac,1)) = 1.0_wrp
        DO i = 2, SIZE(mole_frac, 1)-1
            mole_frac(i) = mole_frac(i-1) + REAL(1.0/interpolate_step_numb, KIND=wrp)
        END DO

        DO i = 1, SIZE(O2C, 1) ! Loop through one compound at a time
            ! Calculate activities
            CALL BAT_properties_calculation_v1(mole_frac, [O2C(i)], [H2C(i)], [Mratio(i)], [BAT_functional_group(i)], &
                & [0.0_wrp], func1, func2, ycalc_water, ycalc_org, activity_water, activity_org, &
                & mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2)
            activity_water(1) = 1.0_wrp ! at x_water = 1.0
            activity_water(interpolate_step_numb + 1_wip) = 0.0_wrp ! at x_water = 0.0
            activity_org(1) = 0.0_wrp ! at x_org = 0.0
            activity_org(interpolate_step_numb + 1_wip) = 1.0_wrp   ! at x_org = 1.0
            
            ! Find water activity at which liquid-liquid separation happens using the calculated activities
            CALL finds_PhaseSep_w_and_org(activity_water, activity_org, phase_sep_check, lower_a_w_sep_index, &
                & upper_a_w_sep_index, matching_Upper_a_w_sep_index)
                        
            IF (ABS(phase_sep_check(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN ! Check if there is phase separation (== 0.0_wrp)
                RH_cross_point(i) = activity_water(upper_a_w_sep_index(1))  ! save phase sep. RH
            ELSE
                RH_cross_point(i) = 0.0_wrp ! no phase separation
            END IF
        END DO

        lower_round = 0.0_wrp
        WHERE (RH_cross_point > 0.0_wrp) lower_round = 1.0_wrp

        not_lower_round = 1.0_wrp
        WHERE (RH_cross_point > 0.0_wrp) not_lower_round = 0.0_wrp

        RH_cross_point = RH_cross_point * lower_round + not_lower_round

        upper_round = 0.0_wrp
        WHERE (RH_cross_point > 1.0_wrp) upper_round = 1.0_wrp

        not_upper_round = 1.0_wrp
        WHERE (RH_cross_point > 1.0_wrp) not_upper_round = 0.0_wrp

        RH_cross_point = not_upper_round * RH_cross_point + upper_round

    END SUBROUTINE biphasic_to_single_phase_RH_master_v4

    !****************************************************************************************
    !*   SUBROUTINE finds_PhaseSep_w_and_org                                                *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Find the indexes of the activity curve that correspond to a region of              *
    !*   liquid-liquid phase separation.                                                    *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE finds_PhaseSep_w_and_org(activity_water, activity_org, phase_sep_check, lower_a_w_sep_index, &
        & upper_a_w_sep_index, matching_Upper_a_w_sep_index)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)     :: activity_water, activity_org ! Water and organic activity
        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)    :: phase_sep_check ! If phase separation exists, == 1.0_wrp, if does not exist = 0.0_wrp
        INTEGER(wip), DIMENSION(:), INTENT(OUT) :: lower_a_w_sep_index, upper_a_w_sep_index, matching_Upper_a_w_sep_index ! Indexes that corresponds to
                                                                                                                          ! the region where liquid-liquid phase separation begins and ends
        ! Local variables
        REAL(wrp), DIMENSION(1)              :: phaseSep_via_activity_curvature_w, phaseSep_via_activity_curvature_org, &
                                                & phaseSep_via_activity_w, phaseSep_via_activity_org
        INTEGER(wip), DIMENSION(1)           :: index_phase_sep_starts_w, index_phase_sep_starts_org, index_phase_sep_end_w, &
                                                & index_phase_sep_end_org
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: activity_water_beta, match_index_prime_a
        REAL(wrp)                            :: match_a_w
        INTEGER(wip), DIMENSION(4)           :: indexes
        INTEGER(wip)                         :: mid_sep_index, m, n, i, match_index_prime_b

        ! Check for phase separation in each activity curve
        ! Organic species
        CALL finds_PhaseSep_and_activity_curve_dips_v2(activity_water, phaseSep_via_activity_w, &
        & phaseSep_via_activity_curvature_w, index_phase_sep_starts_w, index_phase_sep_end_w)
        
        ! Water
        CALL finds_PhaseSep_and_activity_curve_dips_v2(activity_org, phaseSep_via_activity_org, &
        & phaseSep_via_activity_curvature_org, index_phase_sep_starts_org, index_phase_sep_end_org)
        
        ! Start and end indexes of phase separation (for water and organic species)
        indexes = [index_phase_sep_starts_w(1), index_phase_sep_end_w(1), index_phase_sep_starts_org(1), &
            index_phase_sep_end_org(1)]
        IF (ABS(phaseSep_via_activity_curvature_w(1)-1.0_wrp) < tinynumber) THEN
            phase_sep_check = 1.0_wrp

            IF (activity_water(1) < activity_water(SIZE(activity_water, 1))) THEN ! Increasing a_w with index
                lower_a_w_sep_index(1) = MINVAL(indexes)
                upper_a_w_sep_index(1) = MAXVAL(indexes)

                ! Find matching a_w point for upper_a_w_sep_index, as this is not
                ! necessarily the same as lower_a_w_sep_index, and is likely on a metastable
                ! phase-separation branch in Gibbs mix.
                mid_sep_index = FLOOR(0.5_wrp*(lower_a_w_sep_index(1) + upper_a_w_sep_index(1)))
                ALLOCATE(activity_water_beta(mid_sep_index))
                activity_water_beta = activity_water(1:mid_sep_index)
                match_a_w = activity_water(upper_a_w_sep_index(1))

                n = COUNT((activity_water_beta - match_a_w) > 0.0_wrp)
                IF (n > 0) THEN

                    ALLOCATE(match_index_prime_a(n))
                    match_index_prime_a = 0
                    m = 0

                    DO i = 1, mid_sep_index
                        IF ((activity_water_beta(i) - match_a_w) > 0.0_wrp) THEN
                            m = m + 1
                            match_index_prime_a(m) = i
                        END IF
                    END DO
                    matching_Upper_a_w_sep_index(1) = INT(match_index_prime_a(1), KIND = wip)-1;

                ELSE IF (n == 0) THEN
                    match_index_prime_b = MAXLOC(activity_water_beta - match_a_w, 1)
                    matching_Upper_a_w_sep_index(1) = match_index_prime_b - 1
                END IF

            ELSE
                lower_a_w_sep_index(1) = MAXVAL(indexes) ! decreasing a_w with index
                upper_a_w_sep_index(1) = MINVAL(indexes)

                ! Find matching a_w point for upper_a_w_sep_index, as this is not
                ! necessarily the same as lower_a_w_sep_index, and is likely on a metastable
                ! phase-separation branch in Gibbs mix.

                mid_sep_index = FLOOR(0.5_wrp*(lower_a_w_sep_index(1) + upper_a_w_sep_index(1)))

                ALLOCATE(activity_water_beta(SIZE(activity_water, 1) - mid_sep_index + 1))

                activity_water_beta = activity_water(mid_sep_index:SIZE(activity_water, 1))

                match_a_w = activity_water(upper_a_w_sep_index(1))

                n = COUNT(activity_water_beta <= match_a_w)

                ALLOCATE(match_index_prime_a(n))

                m = 0
                DO i = 1, (SIZE(activity_water, 1) - mid_sep_index + 1)
                        IF (activity_water_beta(i) <= match_a_w) THEN
                        m = m + 1
                        match_index_prime_a(m) = i
                    END IF
                END DO

                matching_Upper_a_w_sep_index(1) = mid_sep_index + INT(match_index_prime_a(1), KIND = wip)-1
                END IF
        ELSE
            lower_a_w_sep_index(1) = 1_wip ! no phase sep
            upper_a_w_sep_index(1) = 2_wip
            matching_Upper_a_w_sep_index(1) = 2_wip
            phase_sep_check(1) = 0.0_wrp  ! no phase sep
            
        END IF

    END SUBROUTINE finds_PhaseSep_w_and_org    

    !****************************************************************************************
    !*   SUBROUTINE finds_PhaseSep_and_activity_curve_dips_v2                               *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Builds activity curves to find liquid-liquid phase separation of components.       *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE finds_PhaseSep_and_activity_curve_dips_v2(activity_data1, phaseSep_via_activity, &
        & phaseSep_via_activity_curvature, index_phase_sep_starts, index_phase_sep_end)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)             :: activity_data1                                           ! Activity of organics or water

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)            :: phaseSep_via_activity, phaseSep_via_activity_curvature   ! == 1.0_wrp if phase separation is found
                                                                                                                    ! == 0.0_wrp if phase separation is not found
        INTEGER(wip), DIMENSION(:), INTENT(OUT)         :: index_phase_sep_starts, index_phase_sep_end              ! Indexes to identify composition range that defines
                                                                                                                    ! the miscibility gap
        ! Local variables
        INTEGER(wip)                                    :: L_m, min_index, max_index, index_start, back_index, back_index_temp, activity_data_gap, &
                                                           & restart_match_index, min_index_Idilute, activity_data_gap_start, i
        REAL(wrp), DIMENSION(SIZE(activity_data1, 1)-1) :: activity_data_diff
        REAL(wrp), DIMENSION(SIZE(activity_data1, 1))   :: activity_calc1_diff_sign_change, activity_data_greater_than_1, &
                                                           & activity_calc2_diff_sign_change
        REAL(wrp), DIMENSION(:), ALLOCATABLE            :: not_activity_calc2_diff_sign_change
        REAL(wrp)                                       :: min_value, max_value, min_value_Idilute

        L_m = SIZE(activity_data1, 1)
        DO i = 1, L_m-1
            activity_data_diff(i) = activity_data1(i+1) - activity_data1(i)
		END DO
		
        IF (L_m > 3) THEN

            min_value = MINVAL(activity_data_diff, 1)
            min_index = MINLOC(activity_data_diff, 1)
            max_value = MAXVAL(activity_data_diff, 1)
            max_index = MAXLOC(activity_data_diff, 1)

            ! If the max and min signs are both the same then there is no dip in the activity curve
            IF (ABS(SIGN(1.0_wrp, min_value) - SIGN(1.0_wrp, max_value)) < tinynumber) THEN     ! Try multiplication instead to test sign
                phaseSep_via_activity_curvature(1) = 0.0_wrp
                index_phase_sep_starts(1) = -99999_wip ! NaN
                index_phase_sep_end(1) = -99999_wip ! NaN
            
            ELSE  ! If signs are not the same, there is a dip in the activity curve
                phaseSep_via_activity_curvature(1) = 1.0_wrp    ! Perhaps logical is better
                ! Get mole fraction value at start of phase separation and end of phase separation
                activity_calc1_diff_sign_change = [activity_data_diff(1), activity_data_diff]
                
                ! Find the curve change
                activity_calc2_diff_sign_change = 0.0_wrp
                DO i = 1, L_m
                    IF (ABS(SIGN(1.0_wrp, activity_calc1_diff_sign_change(i)) - SIGN(1.0_wrp, activity_data_diff(1))) < tinynumber) THEN ! If of the same sign
                        activity_calc2_diff_sign_change(i) = 0.0_wrp                   
                    ELSE
                        activity_calc2_diff_sign_change(i) = 1.0_wrp   ! if signs are different 
                    END IF
				END DO 
            
				index_start = L_m
                DO i = 1, L_m
                    IF (activity_calc2_diff_sign_change(i) > 0.0_wrp) THEN
                        index_start = i
                        EXIT
                    END IF
				END DO
				
                ALLOCATE(not_activity_calc2_diff_sign_change(SIZE(activity_calc2_diff_sign_change(index_start:L_m),1)))
                not_activity_calc2_diff_sign_change = activity_calc2_diff_sign_change(index_start:L_m)
				
                DO i = 1, SIZE(not_activity_calc2_diff_sign_change,1)
                    IF (not_activity_calc2_diff_sign_change(i) > 0.0_wrp + tinynumber) THEN
                        not_activity_calc2_diff_sign_change(i) = 0.0_wrp
                    ELSE IF (not_activity_calc2_diff_sign_change(i) < tinynumber) THEN
                        not_activity_calc2_diff_sign_change(i) = 1.0_wrp
                    END IF
                END DO

                DO i = 1, SIZE(not_activity_calc2_diff_sign_change,1)
                    IF (not_activity_calc2_diff_sign_change(i) > 0.0_wrp + tinynumber) THEN
                        back_index_temp = i
                        EXIT
                    END IF
                END DO                
                back_index = index_start - 1_wip + back_index_temp

                IF (back_index < L_m) THEN
                    activity_data_gap = MINLOC(ABS(activity_data1(back_index:L_m)-activity_data1(index_start)), 1)
                    restart_match_index = activity_data_gap + back_index - 1_wip
            
                ELSE
                    restart_match_index = L_m
                END IF
            
                ! Activity change, greater than 1
                activity_data_greater_than_1 = 0.0_wrp
                WHERE (activity_data1 > 1.0_wrp) 
                    activity_data_greater_than_1 = 1.0_wrp
                END WHERE
                
                IF (SUM(activity_data_greater_than_1) > 0.0_wrp) THEN ! greater than 0.0_wrp
                    ! Get min RH in high mole fraction region.
                    min_value_Idilute = MINVAL(activity_data1(index_start:L_m), 1)
                    min_index_Idilute = MINLOC(activity_data1(index_start:L_m), 1)
                    min_index_Idilute = min_index_Idilute + index_start - 1_wip
            
                    ! Get mole fraction in low mole fraction region
                    activity_data_gap_start = MINLOC(ABS(activity_data1(1:index_start)-activity_data1(min_index_Idilute)), 1)
            
                    ! Output
                    ! Check which region starts earlier and ends earlier
                    IF (activity_data_gap_start < index_start) THEN
                        index_phase_sep_starts(1) = activity_data_gap_start
                    ELSE
                        index_phase_sep_starts(1) = index_start
                    END IF
                    IF (min_index_Idilute < restart_match_index) THEN
                        index_phase_sep_end(1) = min_index_Idilute
                    ELSE
                        index_phase_sep_end(1) = restart_match_index
                    END IF
                ELSE
                    !output
                    index_phase_sep_starts(1) = index_start
                    index_phase_sep_end(1) = restart_match_index
                END IF
            END IF
        ELSE ! if L_m < 3
            phaseSep_via_activity = activity_data1
            phaseSep_via_activity_curvature = 0.0_wrp
            index_phase_sep_starts = -99999_wip
            index_phase_sep_end = -99999_wip
        END IF

        activity_data_greater_than_1 = 0.0_wrp
        WHERE (activity_data1 > 1.0_wrp) 
            activity_data_greater_than_1 = 1.0_wrp
        END WHERE
        
        IF (SUM(activity_data_greater_than_1) > 0.0_wrp) THEN ! greater than 0.0_wrp
            phaseSep_via_activity(1) = 1.0_wrp
        
            phaseSep_via_activity(1) = 0.0_wrp
        END IF

    END SUBROUTINE finds_PhaseSep_and_activity_curve_dips_v2

    !****************************************************************************************
    !*   SUBROUTINE single_phase_O2C_point_KGv3                                             *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Sigmoid function that represents the miscibility limit in O/C vs molar mass space. *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE single_phase_O2C_point_KGv3(molarmass_ratio, O2C_single_phase_cross_point)

        IMPLICIT NONE
        
        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: molarmass_ratio
        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: O2C_single_phase_cross_point

        O2C_single_phase_cross_point = 0.205_wrp/(1.0_wrp + EXP(26.6_wrp * (molarmass_ratio - 0.12_wrp)))**0.843_wrp + &
            & 0.225_wrp

    END SUBROUTINE single_phase_O2C_point_KGv3

    !****************************************************************************************
    !*   SUBROUTINE BAT_activity_calc_with_refinement_v1                                    *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Calculates activity organics and water.                                            *
    !*   Calculates the water and organic mass fraction associated with every organic       *
    !*   species (every binary mixture).                                                    *
    !*   This is performed according to the Duhem-Margules relation, which relates the      *
    !*   the molar excess Gibbs energy of mixing to the two mole-fraction-based activity    *
    !*   coefficients (of water and organics). See equations 13 and 14 of BAT paper.        *
    !*                                                                                      *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE BAT_properties_calculation_v1(org_mole_fraction, O2C, H2C, molarmass_ratio, BAT_functional_group, &
        & N2C_values_densityOnly, ln_func1, ln_func2, ycalc1, ycalc2, activity_calc1, activity_calc2, mass_fraction1, &
        & mass_fraction2, Gibbs_RT, dGibbs_RTdx2)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: org_mole_fraction, & ! Mole fraction of organic species
                                                     & O2C, H2C, N2C_values_densityOnly, & ! Elemental ratios of ornanics species (O/C, H/C, N/C)
                                                     & molarmass_ratio ! Molecular weight ratio of organic species (Mwater/Morg)
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: BAT_functional_group ! Main oxygen-bearing functional group of the organic molecule

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: ln_func1, & ! ln(act. coeff. water)
                                                     & ln_func2, & ! ln(act. coeff. org)
                                                     & ycalc1, & ! Mole-fraction based activity coefficient of water
                                                     & ycalc2, & ! Mole-fraction based activity coefficient of organic species
                                                     & activity_calc1, & ! Activity of water
                                                     & activity_calc2, & ! Activity of organic species
                                                     & mass_fraction1, & ! Mass fraction of water
                                                     & mass_fraction2, & ! Mass fraction of organic species
                                                     & Gibbs_RT, & ! Molar excess Gibbs energy of mixing G_excess/RT
                                                     & dGibbs_RTdx2 ! Derivative of molar excess Gibbs energy of mixing G_excess/RT with respect to the mole
                                                                    ! the mole fraction of organic species

        ! Local variables
        TYPE(specialoptions_input)                          :: special_options ! BAT simulation options
        INTEGER(wip)                                        :: allocstat, do_calc, loop_i, n, n1, j
        REAL(wrp)                                           :: tran_lowO2C_fractionOne_phase, tran_lowO2C_sigmoid_bend, tran_lowO2C_sigmoid_shift, &
                                                               & tran_midO2C_sigmoid_bend, tran_midO2C_sigmoid_shift
        REAL(wrp), DIMENSION(SIZE(org_mole_fraction, 1))    :: x2, x2_temp, mole_fraction_mask, not_mole_fraction_mask, &
                                                               & phi2, dphi2dx2, sum1, sum2, dgemix2dx_temp, gemix_temp
        REAL(wrp), DIMENSION(SIZE(O2C, 1))                  :: O2C_temp, molarmass_ratio_temp, Mr_massfrac_final, Mr, densityEst, Onephase_O2C, &
                                                               & mid_transition, O2C_1phase_delta, weight_1phase, weight_2phase, &
                                                               & O2C_1phase_delta_norm, weight_1phase_norm, rhor, scaledMr
        REAL(wrp), DIMENSION(SIZE(org_mole_fraction, 1), 2) :: gemix, dgemix2dx
        CHARACTER(clen)                                     :: fitpar_name1, fitpar_name2
        REAL(wrp), DIMENSION(10)                            :: fitpar_lowO2C, fitpar_midO2C, fitpar_highO2C, fitpar_1phase, fitpar_2phase, fitpar
        REAL(wrp), DIMENSION(2)                             :: coeff

        tran_lowO2C_fractionOne_phase = 0.189974476118418_wrp
        tran_lowO2C_sigmoid_bend = 79.2606902175984_wrp
        tran_lowO2C_sigmoid_shift = 0.0604293454322489_wrp
        tran_midO2C_sigmoid_bend = 75.0159268221068_wrp
        tran_midO2C_sigmoid_shift = 0.000947111285750515_wrp

        ! Shift for equivalent O/C and MW, for activity calculation
        Mr_massfrac_final = molarmass_ratio ! Used in final mass fraction calc and not for activity coefficent calc.
		!PRINT *, O2C, molarmass_ratio, BAT_functional_group
        CALL convert_chemical_structure_to_OH_eqv_v3(O2C, molarmass_ratio, BAT_functional_group, O2C_temp, molarmass_ratio_temp)

        ! Force org mole fraction to be 1 or less
        x2 = org_mole_fraction
        WHERE (org_mole_fraction > 1.0_wrp) x2 = 1.0_wrp
        WHERE (org_mole_fraction < 0.0_wrp) x2 = 0.0_wrp
        
        ! Replace infinite dilution of zero to small value
        WHERE (org_mole_fraction < 1.0E-20_wrp) x2 = 1.0E-20_wrp

        Mr = molarmass_ratio_temp  ! molar1/molar2

        !  Properties
        !  Estimate density of org. components using translated O/C, and molecular weight ratio (to OH-equiv.) and model by Girolami (1994)
        CALL Org_density_Estimate_KGv1((M_water/Mr), O2C_temp, H2C, N2C_values_densityOnly, densityEst)
        ! Estimate O/C limit of miscibility point
        CALL single_phase_O2C_point_KGv3(Mr, Onephase_O2C)

        ! Get region transition properties (BAT model fit parameters for transition from one O/C region to another O/C region).
        mid_transition(1) = Onephase_O2C(1) * 0.75_wrp
        IF (O2C_temp(1) < mid_transition(1)) THEN !  lower to mid O/C region

            ! Data point trasfer weight
            O2C_1phase_delta(1) = O2C_temp(1) - (Onephase_O2C(1) * tran_lowO2C_fractionOne_phase)
            weight_1phase(1) = 1.0_wrp/(1.0_wrp + EXP(-tran_lowO2C_sigmoid_bend * &
            & (O2C_1phase_delta(1) - tran_lowO2C_sigmoid_shift))) ! logistic transfer function 1/(1+e^-(75*x))

            ! Normalize to end point so at mid_transition weight 2 is 1.
            O2C_1phase_delta_norm(1) = O2C_temp(1) - (mid_transition(1) * tran_lowO2C_fractionOne_phase)
            weight_1phase_norm(1) = 1.0_wrp/(1.0_wrp + EXP(-tran_lowO2C_sigmoid_bend * &
                (O2C_1phase_delta_norm(1) - tran_lowO2C_sigmoid_shift))) ! Logistic transfer function 1/(1+e^-(75*x))

            weight_1phase(1) = weight_1phase(1)/weight_1phase_norm(1)
            weight_2phase(1) = 1.0_wrp - weight_1phase(1)

            fitpar_name1 = 'midO2C'
            fitpar_name2 = 'lowO2C'

        ELSE IF (O2C_temp(1) < Onephase_O2C(1) * 2) THEN ! mid to high O/C region

            O2C_1phase_delta(1) = O2C_temp(1) - Onephase_O2C(1)
            weight_1phase(1) = 1.0_wrp/(1.0_wrp + EXP(-tran_midO2C_sigmoid_bend * &
                (O2C_1phase_delta(1) - tran_midO2C_sigmoid_shift)))
            weight_2phase(1) = 1.0_wrp - weight_1phase(1)

            fitpar_name1 = 'highO2C'
            fitpar_name2 = 'midO2C'

        ELSE ! High only region

            fitpar_name1 = 'highO2C'
            fitpar_name2 = 'NaN'

            weight_2phase(1) = 0.0_wrp
            weight_1phase(1) = 1.0_wrp
        
        END IF

        ! Fit properties data
        fitpar_lowO2C = [7.089476E+00_wrp, -0.6226781_wrp, -7.711860E+00_wrp, -1.000000E+02_wrp, -3.885941E+01_wrp, &
            3.081244E-09_wrp, -1.000000E+02_wrp, 6.188812E+01_wrp, -5.988895E+00_wrp, 6.940689E+00_wrp]

        fitpar_midO2C = [5.872214E+00_wrp, -0.9740486_wrp, -4.535007E+00_wrp, -1.000000E+02_wrp, -5.129327E+00_wrp, &
            2.109751E+00_wrp, -2.809232E+01_wrp, -2.367683E+01_wrp, -1.219164E+00_wrp, 4.742729E+00_wrp]

        fitpar_highO2C = [5.921550E+00_wrp, -1.000000E+02_wrp, -2.528295E+00_wrp, -1.000000E+02_wrp, -3.883017E+00_wrp, &
            1.353916E+00_wrp, -7.898128E+00_wrp, -1.160145E+01_wrp, -0.07868187_wrp, 3.650860E+00_wrp]

        ! Get fit parameters in into correct phases: phase 1
        IF (fitpar_name1 == 'highO2C') THEN
            fitpar_1phase = fitpar_highO2C
        ELSE IF (fitpar_name1 == 'midO2C') THEN
            fitpar_1phase = fitpar_midO2C
        ELSE IF (fitpar_name1 == 'lowO2C') THEN
            fitpar_1phase = fitpar_lowO2C
        END IF

        ! Get fit parameters in into correct phases: phase 2
        IF (fitpar_name2 == 'highO2C') THEN
            fitpar_2phase = fitpar_highO2C
        ELSE IF (fitpar_name2 == 'midO2C') THEN
            fitpar_2phase = fitpar_midO2C
        ELSE IF (fitpar_name2 == 'lowO2C') THEN
            fitpar_2phase = fitpar_lowO2C
        END IF

        !! For biphasic line calculations
        !IF (BAT_functional_group(1) == 'initial fitting for biphasic transition') THEN
        !
        !    ! Inital fit to biphasic line
        !    fitpar_1phase = [5.885109E+00_wrp, -9.849013E-01_wrp, -4.731250E+00_wrp, -6.227207E+00_wrp, -5.201652E+00_wrp, &
        !        2.320286E+00_wrp, -3.082297E+01_wrp, -2.584037E+01_wrp, -1.237227E+00_wrp, 4.069905E+00_wrp]
        !    weight_2phase(1) = 0.0_wrp
        !    weight_1phase(1) = 1.0_wrp
        !
        !END IF

        gemix = 0.0_wrp ! Molar excess Gibbs energy of mixing
        dgemix2dx = 0.0_wrp ! Derivative of the molar excess Gibbs energy of mixing with respect to the mole fraction of organic species
        do_calc = 0_wip

        DO loop_i = 1, 2
            IF (loop_i == 1_wip) THEN
                IF (weight_1phase(1) < tinynumber) THEN
                    do_calc = 0_wip
                ELSE
                    fitpar = fitpar_1phase
                    do_calc = 1_wip
                END IF
            ELSE
                IF (weight_2phase(1) < tinynumber) THEN
                    do_calc = 0_wip
                ELSE
                    fitpar = fitpar_2phase
                    do_calc = 1_wip
                END IF
			END IF

            IF (do_calc == 1_wip) THEN
                n = SIZE(fitpar, 1)
                rhor(1) = 0.997_wrp/densityEst(1) ! Assumes water

                ! scaledMr is the scaled molar mass ratio of this mixture's components
                scaledMr(1) = Mr(1) * fitpar(n) * (1.0_wrp + O2C_temp(1))**fitpar(n-1)

                ! phi2 is a scaled volume fraction, which is determined from mole fractions, scaledMr
                ! and an estimate of water/organic density ratio (rhor).
				
                phi2 = x2/(x2 + (1.0_wrp - x2) * scaledMr(1)/rhor(1)) ! and phi1 = 1 - phi2
                dphi2dx2 = (scaledMr(1)/rhor(1)) * (1.0_wrp/(x2 + (1.0_wrp - x2) * scaledMr(1)/rhor(1)))**2 ! the derivative of phi2 with respect to x2 (mole fraction)

                n1 = (n-2)/4
                coeff(1:n1) = fitpar(1:n1) * EXP(fitpar(n1+1:2*n1)*O2C_temp(1)) + fitpar(2*n1+1:3*n1) * &
                    EXP(fitpar(3*n1+1:4*n1)*Mr(1))

                sum1 = 0.0_wrp
                DO j = 1, n1
                    sum1 = sum1 + coeff(j) * (1.0_wrp - 2.0_wrp * phi2)**(j-1)
                END DO

                sum2 = 0.0_wrp
                DO j = 2, n1 ! first derivative
                    sum2 = sum2 + 2.0_wrp * REAL(j-1) * coeff(j) * (1.0_wrp - 2.0_wrp*phi2)**(j-2)
                END DO
            
                gemix(:,loop_i) = phi2 * (1.0_wrp - phi2) * sum1
                dgemix2dx(:,loop_i) = ((1.0_wrp - 2.0_wrp * phi2) * sum1 + phi2 * (1.0_wrp - phi2) * (-sum2)) * dphi2dx2

                END IF
        END DO

        gemix_temp = gemix(:,1) * weight_1phase(1) + gemix(:,2) * weight_2phase(1)
        dgemix2dx_temp = dgemix2dx(:,1) * weight_1phase(1) + dgemix2dx(:,2) * weight_2phase(1)

        ! Calculate the function value funcx1 (= y1(x2)) at point with w2:
        ln_func1 = gemix_temp - x2 * dgemix2dx_temp ! The func value for component 1 = LOG(activity coeff. water)
        ln_func2 = gemix_temp + (1.0_wrp - x2) * dgemix2dx_temp ! The func value of the component 2 = LOG(activity coeff. of the organic)

        WHERE (ln_func1 < -690.7755_wrp) ln_func1 = -690.7755_wrp
        WHERE (ln_func1 > 690.7755_wrp) ln_func1 = 690.7755_wrp
        ycalc1 = EXP(ln_func1) ! activity coefficient water

        WHERE (ln_func2 < -690.7755_wrp) ln_func2 = -690.7755_wrp
        WHERE (ln_func2 > 690.7755_wrp) ln_func2 = 690.7755_wrp
        ycalc2 = EXP(ln_func2) ! activity coefficient org

        activity_calc1 = ycalc1 * (1.0_wrp - x2) ! activity water
        activity_calc2 = ycalc2 * x2 ! avtivity org

        mass_fraction1 = (1.0_wrp - x2) * Mr_massfrac_final(1)/((1.0_wrp - x2) * (Mr_massfrac_final(1) - 1.0_wrp) + 1.0_wrp) ! mass fraction water
        mass_fraction2 = 1.0_wrp - mass_fraction1 ! mass fraction org
        Gibbs_RT = gemix_temp ! G_excess/RT
        dGibbs_RTdx2 = dgemix2dx_temp ! d(G_excess/RT)/d (mole fraction org)

    END SUBROUTINE BAT_properties_calculation_v1

    !****************************************************************************************
    !*   SUBROUTINE BAT_activity_calc_with_refinement_v1                                    *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Calculates activity of organics and water.                                         *
    !*   Calculates the water and organic mass fraction associated with every organic       *
    !*   species (every binary mixture).                                                    *
    !*   A refinement via 501 mole fraction interpolation points is done when RH is greater *
    !*   than 90%. The mole fraction of organic species above 90 % relative humidity        *
    !*   is determined through the refinement method so that water activity in the gas      *
    !*   phase (RH) matches the water activity in the liquid phase.                         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE BAT_activity_calc_with_refinement_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio, &
        & BAT_functional_group, refinement_mode, aw_desired, N2C_values_densityOnly, func1, func2, ycalc_water, &
        & ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2, &
        & mole_frac_fit, error_out)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: O2C_values, H2C_values, N2C_values_densityOnly, & ! Elemental ratios (O/C, H/C, N/C) of org.
                                                     & molarmass_ratio, & ! Molecular weight ratio of org. (Mwater/Morg)
                                                     & mole_frac, & ! Mole fraction of org.
                                                     & aw_desired ! Water activity the model is evaluated at
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: BAT_functional_group, & ! Main oxygen-bearing functional group of the organic molecule
                                                     & refinement_mode ! Use interpolation or not

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: mass_fraction1, mass_fraction2, & ! Mass fraction of water (=1) and org. (=2)
                                                     & Gibbs_RT, dGibbs_RTdx2, & ! Molar excess Gibbs energy (G/RT) and its derivative w.r.t. org. mole frac.
                                                     & func1, func2, & ! ln(act. coeff. water) and ln(act. coeff. org.)
                                                     & ycalc_water, ycalc_org, & ! act. coeff. water and act. coeff. org.
                                                     & activity_water, activity_org, &! activity water and activity org.
                                                     & error_out, mole_frac_fit ! relative error in water activity: target vs calculated; corresponding org.
                                                                                ! mole fraction at equilibrium

        ! Local variables
        INTEGER(wip)               :: i, aw_index, interpolate_step_numb
        INTEGER(wip), DIMENSION(1) :: index_phase_sep_starts, index_phase_sep_end
        REAL(wrp), DIMENSION(1)    :: phaseSep_via_activity, phaseSep_via_activity_curvature
        REAL(wrp), DIMENSION(501)  :: mole_frac_interpolation, func1_interpolation, func2_interpolation, &
                                      & ycal_water_interpolation, ycalc_org_interpolation, activity_water_interpolation, &
                                      & activity_org_interpolation, mass_fraction1_interpolation, mass_fraction2_interpolation, &
                                      & Gibbs_RT_interpolation, dGibbs_RTdx2_interpolation

        interpolate_step_numb = 500_wip ! Interpolate number of steps the resolution interpolate uses
        mole_frac_fit = 0.0_wrp

        IF (refinement_mode(1) == 'none') THEN
            CALL BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio,  &
            & BAT_functional_group, N2C_values_densityOnly, func1, func2,  &
            & ycalc_water, ycalc_org,  activity_water, activity_org, mass_fraction1, mass_fraction2, Gibbs_RT,  dGibbs_RTdx2)

            error_out = -99999.0_wrp
            mole_frac_fit = mole_frac

        ! Interpolate using the mole fraction in phase beta (organic-rich phase)
        ELSE IF (refinement_mode(1) == 'interpolatebeta') THEN

            mole_frac_interpolation = 0.0_wrp
            DO i = 2, SIZE(mole_frac_interpolation, 1)
                mole_frac_interpolation(i) = mole_frac_interpolation(i-1) + REAL(1.0/interpolate_step_numb, KIND=wrp)
            END DO
            mole_frac_interpolation(SIZE(mole_frac_interpolation, 1)) = 1.0_wrp

            ! Calculate a water activity curve using 501 values of mole fraction
            CALL BAT_properties_calculation_v1(mole_frac_interpolation, O2C_values, H2C_values, molarmass_ratio, &
                & BAT_functional_group, N2C_values_densityOnly, func1_interpolation, func2_interpolation, &
                & ycal_water_interpolation, ycalc_org_interpolation, activity_water_interpolation, activity_org_interpolation, &
                & mass_fraction1_interpolation, mass_fraction2_interpolation, Gibbs_RT_interpolation, dGibbs_RTdx2_interpolation)

            ! Check for liquid-liquid phase separation in water activity curve
            CALL finds_PhaseSep_and_activity_curve_dips_v2(activity_water_interpolation, phaseSep_via_activity, &
            & phaseSep_via_activity_curvature, index_phase_sep_starts, index_phase_sep_end)

            ! If phase separation is found
            IF (ABS(phaseSep_via_activity_curvature(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN
                IF (index_phase_sep_end(1) < SIZE(activity_water_interpolation, 1)) THEN
                    ! Apply linear interpolation if the last mole fraction of phase separation region is not equal to the last mole fraction of the curve
                    CALL interp1(activity_water_interpolation(index_phase_sep_end(1):SIZE(activity_water_interpolation, 1)), &
                    & mole_frac_interpolation(index_phase_sep_end(1):SIZE(activity_water_interpolation, 1)), aw_desired, mole_frac_fit)
                ELSE
                    mole_frac_fit = mole_frac_interpolation(SIZE(mole_frac_interpolation, 1))
                END IF
            ELSE ! If phase separation not found
                ! Apply linear interpolation using entire aw curve
                CALL interp1(activity_water_interpolation(1:SIZE(activity_water_interpolation, 1)), &
                & mole_frac_interpolation(1:SIZE(activity_water_interpolation, 1)), aw_desired, mole_frac_fit)
            END IF

            IF (mole_frac_fit(1) < -90000.0_wrp) THEN ! == -99999.0_wrp
                mole_frac_fit(1) = mole_frac_interpolation(index_phase_sep_end(1))
            END IF

            aw_index = MINLOC(ABS(mole_frac_interpolation - mole_frac_fit(1)), 1)

            !  Final output
            CALL BAT_properties_calculation_v1(mole_frac_fit, O2C_values, H2C_values, molarmass_ratio,  &
            & BAT_functional_group, N2C_values_densityOnly, func1, func2,  ycalc_water, &
            & ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2)

            ! Absolute error in water activity after refinement: aw target vs aw calculated
            error_out = ABS(aw_desired-activity_water) !ABS(aw_desired-activity_water)/aw_desired
			
        ! Interpolate using the mole fraction in phase alpha
        ELSE IF (refinement_mode(1) == 'interpolatealpha') THEN

            mole_frac_interpolation = 0.0_wrp
            DO i = 2, SIZE(mole_frac_interpolation, 1)
                mole_frac_interpolation(i) = mole_frac_interpolation(i-1) + REAL(1.0/interpolate_step_numb, KIND=wrp)
            END DO
            mole_frac_interpolation(SIZE(mole_frac_interpolation,1)) = 1.0_wrp
            
            ! Calculate a water activity curve using 501 values of mole fraction
            CALL BAT_properties_calculation_v1(mole_frac_interpolation, O2C_values, H2C_values, molarmass_ratio, &
            & BAT_functional_group, N2C_values_densityOnly, func1_interpolation, func2_interpolation, &
            & ycal_water_interpolation, ycalc_org_interpolation, activity_water_interpolation, activity_org_interpolation, &
            & mass_fraction1_interpolation, mass_fraction2_interpolation, Gibbs_RT_interpolation, dGibbs_RTdx2_interpolation)
            
            ! Check for liquid-liquid phase separation in water activity curve
            CALL finds_PhaseSep_and_activity_curve_dips_v2(activity_water_interpolation, phaseSep_via_activity, &
            & phaseSep_via_activity_curvature, index_phase_sep_starts, index_phase_sep_end)

            ! If phase separation is found
            IF (ABS(phaseSep_via_activity_curvature(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN
                IF (index_phase_sep_end(1) < SIZE(activity_water_interpolation, 1)) THEN
                    ! Apply linear interpolation if the last mole fraction of phase separation region is not equal to the last mole fraction of the curve
                    CALL interp1(activity_water_interpolation(1:index_phase_sep_starts(1)), &
                    & mole_frac_interpolation(1:index_phase_sep_starts(1)), aw_desired, mole_frac_fit)
                ELSE
                    mole_frac_fit = mole_frac_interpolation(SIZE(mole_frac_interpolation, 1))
                END IF
            ELSE ! If no phase separation is found
                ! Apply linear interpolation using the entire water activity curve
                CALL interp1(activity_water_interpolation(1:SIZE(mole_frac_interpolation, 1)), &
                & mole_frac_interpolation(1:SIZE(mole_frac_interpolation, 1)), aw_desired, mole_frac_fit)
            END IF

            IF (mole_frac_fit(1) < -90000.0_wrp) THEN ! == -99999.0_wrp
                mole_frac_fit(1) = mole_frac_interpolation(index_phase_sep_starts(1))
            END IF

            aw_index = MINLOC(ABS(mole_frac_interpolation - mole_frac_fit(1)), 1)

            ! Final output
            CALL BAT_properties_calculation_v1(mole_frac_fit, O2C_values, H2C_values, molarmass_ratio, &
            & BAT_functional_group, N2C_values_densityOnly, func1, func2, ycalc_water, &
            & ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2)
            
            ! Absolute error in water activity after refinement: aw target vs aw calculated
            error_out = ABS(aw_desired-activity_water) !ABS(aw_desired-activity_water)/aw_desired

        ELSE
            PRINT *, 'pick BAT refinement method'
        END IF
		
    END SUBROUTINE BAT_activity_calc_with_refinement_v1

    !****************************************************************************************
    !*   SUBROUTINE q_alpha_transfer_vs_aw_calc_v1                                          *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Makes a squeezed logistic function to transfer from q_alpha ~0 to q_alpha ~1.      *
    !*   Uses a sigmoidal function to approximate the fractional partitioning of a          *
    !*   component to phase alpha (water-rich phase).                                       *
    !*   (Equation 21 of BAT paper)                                                         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE q_alpha_transfer_vs_aw_calc_v1(a_w_sep, aw_series, VBSBAT_options, q_alpha_value)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:,:), INTENT(IN)  :: a_w_sep, aw_series ! water activity at which liquid-liquid phase separation occurs and
                                                                     ! input water activity (aw which the model is evaluated)
        TYPE(VBSBAT_options_input), INTENT(IN) :: VBSBAT_options     ! Simulation options

        ! Outputs
        REAL(wrp), DIMENSION(:,:), INTENT(OUT) :: q_alpha_value      ! Fraction of each organic species in liquid phase alpha (w.r.t total j in
                                                                     ! liquid phases
        ! Local variables
        REAL(wrp), DIMENSION(SIZE(a_w_sep,1),SIZE(a_w_sep, 2)) :: mask_of_miscible_points, not_mask_of_miscible_points, &
                                                                  & delta_a_w_sep, above_min_delta_a_w_sep_value, &
                                                                  & not_above_min_delta_a_w_sep_value, sigmoid_curve_parameter, q_alpha_value_temp
        INTEGER(wip)                                           :: i,j

        ! Values held for correction at the end
        mask_of_miscible_points = 0.0_wrp
        not_mask_of_miscible_points = 0.0_wrp
        WHERE (a_w_sep < 0.0 + tinynumber) mask_of_miscible_points = 1.0_wrp ! == 0
        WHERE (mask_of_miscible_points < 0.0 + tinynumber) not_mask_of_miscible_points = 1.0_wrp

        ! Spread (aw gap from aw at liquid-liquid separation point to complete aqueous dilution) in transfer from 50/50 point (sigmoid half-width)
        delta_a_w_sep = 1.0_wrp - a_w_sep

        ! Check min value allowed
        above_min_delta_a_w_sep_value = delta_a_w_sep*0.0_wrp
        not_above_min_delta_a_w_sep_value = delta_a_w_sep*0.0_wrp
        WHERE (delta_a_w_sep > VBSBAT_options%q_alpha%min_spread_in_aw) above_min_delta_a_w_sep_value = 1.0_wrp
        WHERE (above_min_delta_a_w_sep_value < 0.0 + tinynumber) not_above_min_delta_a_w_sep_value = 1.0_wrp
            
        delta_a_w_sep = delta_a_w_sep * above_min_delta_a_w_sep_value + not_above_min_delta_a_w_sep_value * &
            & VBSBAT_options%q_alpha%min_spread_in_aw

        ! Calculate curve parameter of sigmoid
        sigmoid_curve_parameter = log((1.0_wrp/(1.0_wrp-VBSBAT_options%q_alpha%q_alpha_at_1phase_aw) - &
            1.0_wrp))/delta_a_w_sep


        IF(ABS((1.0_wrp/(1.0_wrp-VBSBAT_options%q_alpha%q_alpha_at_1phase_aw)-1.0_wrp)) < tinynumber) THEN
            sigmoid_curve_parameter = (1.0_wrp/(1.0_wrp-VBSBAT_options%q_alpha%q_alpha_at_1phase_aw)-1.0_wrp)/delta_a_w_sep
        ELSE
            sigmoid_curve_parameter = log((1.0_wrp/(1.0_wrp-VBSBAT_options%q_alpha%q_alpha_at_1phase_aw) - &
                1.0_wrp))/delta_a_w_sep
        END IF

        ! Calculate q_alpha value
        q_alpha_value_temp = (sigmoid_curve_parameter * (aw_series - a_w_sep + delta_a_w_sep))
        WHERE(q_alpha_value_temp < -690.7755_wrp)
            q_alpha_value_temp = -690.7755_wrp
        ELSEWHERE (q_alpha_value_temp > 690.7755_wrp)
            q_alpha_value_temp = 690.7755_wrp
        END WHERE

        ! q_alpha of organic species as a function of water activity
        q_alpha_value = 1.0_wrp - 1.0_wrp/(1.0_wrp + exp(sigmoid_curve_parameter * (aw_series - a_w_sep + delta_a_w_sep)))

        ! Apply mask for complete miscibility, turns miscible organics to q_alpha=1 for all water activities a_w
        q_alpha_value = q_alpha_value * not_mask_of_miscible_points + mask_of_miscible_points

    END SUBROUTINE q_alpha_transfer_vs_aw_calc_v1

END MODULE BAT_MOD
