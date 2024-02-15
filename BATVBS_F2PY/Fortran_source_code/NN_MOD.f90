!****************************************************************************************
!*   MODULE NN_MOD                                                                      *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   This module contains two main subroutines:                                         *
!*   1) BAT Neural Network, which estimates the mole fraction of organics at a given RH.*
!*   2) VBS Neural Network, which estimates the partitioning coefficients of organics.  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 06-05-2021                                                      *
!****************************************************************************************

MODULE NN_MOD 

    USE PRECISION_MOD ! Module that defines the precision/length of REAL/INTERGER/CHAR data type
    USE DATATYPE_MOD ! Module that defines derived data types
    USE BAT_MOD ! Module that contains BAT subroutines
    USE TOOLS_MOD ! Module that contains subroutines to perform different recurrent calculations
    USE NN_DATA_MOD ! Module that contains all the constants used in the neural network simulations

    IMPLICIT NONE

CONTAINS

    !****************************************************************************************
    !*   SUBROUTINE inverted_NNBAT_v8                                                       *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that estimates the mole fraction of organics at a given relative        *
    !*   humidity (aw) using their elemental ratios (O/C and H/C) and molecular weight      *
    !*   (via the BAT-Neural Network).                                                      *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************
        
    SUBROUTINE inverted_NNBAT_v8(O2C, H2C, Mratio, a_w, BAT_functional_group, mole_frac_org_alpha, mole_frac_org_beta)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)                :: O2C, H2C, Mratio, a_w ! Elemental O/C and H2C ratios,
                                                                                    ! molecular weight ratio (Mwater/Morg), water activity
        CHARACTER(clen), DIMENSION(:), INTENT(IN)          :: BAT_functional_group  ! Type of oxygen-bearing functional group

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)               :: mole_frac_org_alpha, mole_frac_org_beta ! mole frac. of organics in phases
                                                                                                      ! alpha and beta
        ! Local variables
        REAL(wrp), DIMENSION(SIZE(O2C, 1))                 :: O2C_temp, Mratio_temp
        REAL(wrp), DIMENSION(SIZE(mole_frac_org_alpha, 1)) :: lower_round, upper_round, not_lower_round, not_upper_round
        REAL(wrp), DIMENSION(1)                            :: mole_frac_org_beta1, mole_frac_org_beta2, mole_frac_org_alpha1, &
                                                              & mole_frac_org_alpha2, O2C_single_phase_cross_point, weight_lower, &
                                                              & weight_higher, transfert, RH_cross_point, mole_frac_org_beta_i, &
                                                              & mole_frac_org_alpha_i
        REAL(wrp), DIMENSION(3)                            :: x1, x_alpha
        INTEGER(wip)                                       :: i

        ! Translate the properties of functionalized molecules to a hypothetical OH-equivalent molecule of modified O/C and molecular weight
        CALL convert_chemical_structure_to_OH_eqv_v3(O2C, Mratio, BAT_functional_group, O2C_temp, Mratio_temp)

        mole_frac_org_alpha = 0.0_wrp
        mole_frac_org_beta = 0.0_wrp

        DO i = 1, SIZE(O2C, 1)

            ! Predict the limit of miscibility
            CALL single_phase_O2C_point_KGv3([Mratio_temp(i)], O2C_single_phase_cross_point)
            x1(1) = O2C_temp(i)
            x1(2) = Mratio_temp(i)
            x1(3) = a_w(i)
            CALL invert_NN_85to90_transfer_weights([a_w(i)], transfert, weight_lower, weight_higher)
            
            IF (O2C_temp(i) < O2C_single_phase_cross_point(1)) THEN ! biphasic

                CALL biphasic_to_single_phase_RH_master_v4([O2C_temp(i)], [H2C(i)], [Mratio_temp(i)], &
                    & [BAT_functional_group(i)], RH_cross_point)
                
                IF (ABS(weight_lower(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_biphasic_orgPhase_KGv4(x1, mole_frac_org_beta_i) ! org phase
                    mole_frac_org_beta(i) = mole_frac_org_beta_i(1)
                ELSE IF (ABS(weight_higher(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN
     
                    CALL inverted_NNBAT_biphasic_orgPhase85to1_KGv4(x1, mole_frac_org_beta_i) ! org phase
                    mole_frac_org_beta(i) = mole_frac_org_beta_i(1)
                ELSE IF (ABS(transfert(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_biphasic_orgPhase_KGv4(x1, mole_frac_org_beta1)
                    CALL inverted_NNBAT_biphasic_orgPhase85to1_KGv4(x1, mole_frac_org_beta2)
                    mole_frac_org_beta(i) = weight_lower(1) * mole_frac_org_beta1(1) + weight_higher(1) * &
                        mole_frac_org_beta2(1)
                ELSE

                    PRINT *, 'invert error'
                    STOP
                END IF

                IF (a_w(i) > RH_cross_point(1)) THEN ! last_RH_cross_point

                    x_alpha = x1
                ELSE

                    x_alpha(1) = O2C_temp(i)
                    x_alpha(2) = Mratio_temp(i)
                    x_alpha(3) = RH_cross_point(1)
                END IF
                
                IF (ABS(weight_lower(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_biphasic_aqPhase_KGv4(x_alpha, mole_frac_org_alpha_i)
                    mole_frac_org_alpha(i) = mole_frac_org_alpha_i(1)
                ELSE IF (ABS(weight_higher(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN
 
                    CALL inverted_NNBAT_biphasic_aqPhase85to1_KGv4(x_alpha, mole_frac_org_alpha_i)
                    mole_frac_org_alpha(i) = mole_frac_org_alpha_i(1) 
                ELSE IF (ABS(transfert(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_biphasic_aqPhase_KGv4(x_alpha, mole_frac_org_alpha1)
                    CALL inverted_NNBAT_biphasic_aqPhase85to1_KGv4(x_alpha, mole_frac_org_alpha2)
                    mole_frac_org_alpha(i) = weight_lower(1) * mole_frac_org_alpha1(1) + weight_higher(1) * &
                        mole_frac_org_alpha2(1)

                ELSE
    
                    PRINT *, 'invert error'
                    STOP
                END IF

            ELSE ! single phase

                IF (ABS(weight_lower(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_singlephase_only_KGv4(x1, mole_frac_org_alpha_i)
                    mole_frac_org_alpha(i) = mole_frac_org_alpha_i(1)
                    
                ELSE IF (ABS(weight_higher(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_singlephase_only85to1_KGv4(x1, mole_frac_org_alpha_i)
                    mole_frac_org_alpha(i) = mole_frac_org_alpha_i(1)

                ELSE IF (ABS(transfert(1)-REAL(1.0, KIND = wrp)) < tinynumber) THEN

                    CALL inverted_NNBAT_singlephase_only_KGv4(x1, mole_frac_org_alpha1)
                    CALL inverted_NNBAT_singlephase_only85to1_KGv4(x1, mole_frac_org_alpha2)
                    mole_frac_org_alpha(i) = weight_lower(1) * mole_frac_org_alpha1(1) + weight_higher(1) * &
                        & mole_frac_org_alpha2(1)
                    
                ELSE

                    PRINT *, 'invert error'
                    STOP
                END IF

                mole_frac_org_beta(i) = mole_frac_org_alpha(i)

            END IF

        END DO
        
        ! round to zero alpha
        lower_round = 0.0_wrp
        WHERE (mole_frac_org_alpha > 0.0_wrp) lower_round = 1.0_wrp

        not_lower_round = 1.0_wrp
        WHERE (mole_frac_org_alpha > 0.0_wrp) not_lower_round = 0.0_wrp

        mole_frac_org_alpha = mole_frac_org_alpha * lower_round + not_lower_round * 10.0_wrp**(-8)

        ! round max to 1
        upper_round = 0.0_wrp
        WHERE(mole_frac_org_alpha > 1.0_wrp) upper_round = 1.0_wrp

        not_upper_round = 1.0_wrp
        WHERE(mole_frac_org_alpha > 1.0_wrp) not_upper_round = 0.0_wrp

        mole_frac_org_alpha = not_upper_round * mole_frac_org_alpha + upper_round

        ! round to zero
        lower_round = 0.0_wrp
        WHERE (mole_frac_org_beta > 0.0_wrp) lower_round = 1.0_wrp

        not_lower_round = 1.0_wrp
        WHERE (mole_frac_org_beta > 0.0_wrp) not_lower_round = 0.0_wrp

        mole_frac_org_beta = mole_frac_org_beta * lower_round + not_lower_round * 10.0_wrp**(-8)

        ! round max to 1
        upper_round = 0.0_wrp
        WHERE(mole_frac_org_beta > 1.0_wrp) upper_round = 1.0_wrp

        not_upper_round = 1.0_wrp
        WHERE(mole_frac_org_beta > 1.0_wrp) not_upper_round = 0.0_wrp

        mole_frac_org_beta = not_upper_round * mole_frac_org_beta + upper_round

    END SUBROUTINE inverted_NNBAT_v8

    !****************************************************************************************
    !*   SUBROUTINE VBSBAT_neural_network_estimate_v1                                       *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that estimates the partitioning coefficient Ej of each organic species  *
    !*   in the system (via the VBS-Neural Network).                                        *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE VBSBAT_neural_network_estimate_v1(Cj_sat_ugPm3, Cj_ugPm3, O2C, molecular_weight_gPmol, &
        & mass_fraction_water_beta, a_water, VBSBAT_NN_options, partition_coefficent_estimate)

        IMPLICIT NONE

        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)     :: Cj_sat_ugPm3, & ! Saturation mass concentration of each org.
                                                   & Cj_ugPm3 , O2C, molecular_weight_gPmol, & ! Mass concentration, O/C ratio and molecular weight of each org.
                                                   & mass_fraction_water_beta, a_water ! Mass fraction of water in phase beta (per org.) and water activity (RH)
        TYPE(VBSBAT_NN_options_input)           :: VBSBAT_NN_options ! Different options for the neural network.

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT)    :: partition_coefficent_estimate ! Ej guess
                                                                                 ! Neural network guess of the fraction of each org. j
                                                                                 ! that partitioned to the condensed phase
        ! Local variables
        INTEGER(wip)                            :: Sc, min_val, max_val, i, j, p, input_i
        INTEGER(wip), DIMENSION(:), ALLOCATABLE :: index_to_bin, unique, bins_to_use
        REAL(wrp), DIMENSION(:), ALLOCATABLE    :: Csat_j_bins, input_Cj, input_O2C, input_M, input_mf_water, &
                                                   & input, pc_temp, O2C_miscibility_line, O2C_miscibility_delta

        ! Blank output
        Sc = SIZE(Cj_ugPm3, 1)
        partition_coefficent_estimate = 0.0_wrp

        ! Bining of inputs
        ALLOCATE(Csat_j_bins(12))
        Csat_j_bins = [1.0E+03_wrp * 0.0000_wrp, 1.0E+03_wrp * 0.000000003162278_wrp, 1.0E+03_wrp * 0.000000031622777_wrp, &
            1.0E+03_wrp * 0.000000316227766_wrp, 1.0E+03_wrp * 0.000003162277660_wrp, 1.0E+03_wrp * 0.000031622776602_wrp, &
            1.0E+03_wrp * 0.000316227766017_wrp, 1.0E+03_wrp * 0.003162277660168_wrp, 1.0E+03_wrp * 0.031622776601684_wrp, &
            1.0E+03_wrp * 0.316227766016838_wrp, 1.0E+03_wrp * 3.162277660168380_wrp , HUGE(999999999.99_wrp)] ! shift to bin
        ! edges and extend to full numberline 0 to inf

        ! Bin indexes
        ALLOCATE(index_to_bin(Sc))

        DO i=1, SIZE(Cj_sat_ugPm3, 1)
            DO j=1, SIZE(Csat_j_bins, 1) - 1 ! there are 11 bins
                IF (j < SIZE(Csat_j_bins, 1) - 1) THEN
                    IF (Csat_j_bins(j) <= Cj_sat_ugPm3(i) .AND. Cj_sat_ugPm3(i) < Csat_j_bins(j+1)) THEN
                        index_to_bin(i) = j
                    END IF
                ELSE IF (j == SIZE(Csat_j_bins, 1) - 1) THEN  ! for the bin # 11
                    IF (Csat_j_bins(j) <= Cj_sat_ugPm3(i) .AND. Cj_sat_ugPm3(i) <= Csat_j_bins(j+1)) THEN
                        index_to_bin(i) = j
                    END IF
                END IF
            END DO
        END DO

        IF (SIZE(O2C,1) > 1) THEN
            ! find unique
            ALLOCATE(unique(SIZE(index_to_bin,1)))

            min_val = MINVAL(index_to_bin)-1
            max_val = MAXVAL(index_to_bin)
            
            i = 0
            DO WHILE (min_val < max_val)
                i = i + 1
                min_val = MINVAL(index_to_bin, MASK = index_to_bin > min_val)
                unique(i) = min_val
            END DO
            
            ALLOCATE(bins_to_use(i), source=unique(1:i))

        ELSEIF (SIZE(O2C,1) == 1) THEN
            ALLOCATE(bins_to_use(1))
            bins_to_use = SUM(index_to_bin,1)
        END IF

        ALLOCATE(input_Cj(SIZE(Csat_j_bins,1)-1), input_O2C(SIZE(Csat_j_bins,1)-1), input_M(SIZE(Csat_j_bins,1)-1), &
            input_mf_water(SIZE(Csat_j_bins,1)-1))

        input_Cj = 0.0_wrp
        input_O2C = 0.0_wrp
        input_M = 0.0_wrp
        input_mf_water = 0.0_wrp

        ! sort data into bins
        IF (SIZE(bins_to_use,1) > 1) THEN
            DO i = 1, SIZE(bins_to_use,1)
                input_i = bins_to_use(i)
                input_Cj(input_i) = SUM(Cj_ugPm3, DIM = 1, MASK = index_to_bin == input_i)
                input_O2C(input_i) = SUM(O2C, DIM = 1, MASK = index_to_bin == input_i)/(COUNT(index_to_bin == input_i))
                input_M(input_i) = SUM(molecular_weight_gPmol,DIM = 1, MASK = index_to_bin == input_i) / &
                    (COUNT(index_to_bin == input_i))
                input_mf_water(input_i) = SUM(mass_fraction_water_beta,DIM = 1, MASK = index_to_bin == input_i) / &
                    (COUNT(index_to_bin == input_i))
            END DO
            
        ELSE IF (SIZE(bins_to_use,1) == 1) THEN
            input_Cj = Cj_ugPm3
            input_O2C = O2C
            input_M = molecular_weight_gPmol
            input_mf_water = mass_fraction_water_beta
        END IF


        IF (VBSBAT_NN_options%NN_type == 'mean_properties') THEN

            ! could of only used 1 input for water mass fraction, may update NN later
            ALLOCATE(input(SIZE(input_Cj,1) + SIZE(input_mf_water,1) + SIZE(a_water,1) + 2))
            input = [input_Cj, SUM(input_O2C,1)/SIZE(input_O2C,1), SUM(input_M,1)/SIZE(input_M,1), input_mf_water, a_water]

            ALLOCATE(pc_temp(11))
            CALL NN_VBSBAT_layers24_meanProp_v1(input, pc_temp)

        ELSE IF (VBSBAT_NN_options%NN_type == 'individual_properties') THEN

            ALLOCATE(input(SIZE(input_Cj,1) + SIZE(input_O2C,1) + SIZE(input_M,1) + SIZE(input_mf_water,1) + SIZE(a_water,1)))
            input = [input_Cj, input_O2C, input_M, input_mf_water, a_water]

            ALLOCATE(pc_temp(11))
            CALL NN_VBSBAT_layers28_individProp_v1(input, pc_temp)

        ELSE IF (VBSBAT_NN_options%NN_type == 'individual_properties_v2') THEN
            ALLOCATE(O2C_miscibility_line(SIZE(input_M, 1)))

            CALL single_phase_O2C_point_KGv3(M_water/input_M, O2C_miscibility_line)
            O2C_miscibility_delta = input_O2C - O2C_miscibility_line

            WHERE (input_M < tinynumber) input_M = 0.5_wrp * M_water
            WHERE (input_Cj > 0.0_wrp + tinynumber .OR. input_Cj < 0.0_wrp - tinynumber) input_Cj = log10(input_Cj)


            ALLOCATE(input( SIZE(input_Cj,1)+SIZE(O2C_miscibility_delta,1)+SIZE(input_M,1) + &
                SIZE(input_mf_water,1)+SIZE(a_water,1)))

            input = [input_Cj, O2C_miscibility_delta, M_water/input_M, input_mf_water, a_water]

            ALLOCATE(pc_temp(11))
            CALL NN_VBSBAT_layers20_individProp_v2(input, pc_temp)
        ELSE
            PRINT *, 'Select VBSBAT_NN_options%NN_type'
        END IF

        WHERE (pc_temp < 0.0_wrp)
            pc_temp = 0.0_wrp
        ELSEWHERE (pc_temp > 1.0_wrp)
            pc_temp = 1.0_wrp
        END WHERE

        ! Resort to match inputs
        partition_coefficent_estimate =  Cj_sat_ugPm3 * 0.0_wrp

        ! Sort data into bins
        DO i = 1, SIZE(bins_to_use, 1)
            input_i = bins_to_use(i)
            DO j = 1, SIZE(index_to_bin, 1)
                IF (index_to_bin(j) == input_i) THEN
                    partition_coefficent_estimate(j) = pc_temp(input_i)
                END IF
            END DO
        END DO

    END SUBROUTINE VBSBAT_neural_network_estimate_v1

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*    All the following subroutines are related to the layers of artificial neurons     *
!*    used in the deep belief network.                                                  *
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

    SUBROUTINE inverted_NNBAT_singlephase_only85to1_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(30)             :: b1
        REAL(wrp), DIMENSION(90)             :: IW1_1_temp
        REAL(wrp), DIMENSION(30,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,30)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_singlephase_only85to1
        x1_step1%gain = x1_step1_gain_singlephase_only85to1
        x1_step1%ymin = x1_step1_ymin_singlephase_only85to1

        ! Layer 1
        b1 = b1_singlephase_only85to1
        IW1_1_temp = IW1_1_temp_singlephase_only85to1

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 =  b2_singlephase_only85to1
        LW2_1 = LW2_1_singlephase_only85to1


        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_singlephase_only85to1
        y1_step1%gain =  y1_step1_gain_singlephase_only85to1
        y1_step1%xoffset = y1_step1_xoffset_singlephase_only85to1

        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1 ! samples

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_singlephase_only85to1_KGv4

!****************************************************************************************

    SUBROUTINE inverted_NNBAT_singlephase_only_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(20)             :: b1
        REAL(wrp), DIMENSION(60)             :: IW1_1_temp
        REAL(wrp), DIMENSION(20,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,20)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip) :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_singlephase_only
        x1_step1%gain =  x1_step1_gain_singlephase_only
        x1_step1%ymin = x1_step1_ymin_singlephase_only

        ! Layer 1
        b1 = b1_singlephase_only
        IW1_1_temp = IW1_1_temp_singlephase_only
        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 =  b2_singlephase_only
        LW2_1 = LW2_1_singlephase_only

        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_singlephase_only
        y1_step1%gain = y1_step1_gain_singlephase_only
        y1_step1%xoffset = y1_step1_xoffset_singlephase_only


        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1 ! samples

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_singlephase_only_KGv4

!****************************************************************************************

    SUBROUTINE inverted_NNBAT_biphasic_aqPhase85to1_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(30)             :: b1
        REAL(wrp), DIMENSION(90)             :: IW1_1_temp
        REAL(wrp), DIMENSION(30,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,30)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_biphasic_aqPhase85to1
        x1_step1%gain = x1_step1_gain_biphasic_aqPhase85to1
        x1_step1%ymin = x1_step1_ymin_biphasic_aqPhase85to1

        ! Layer 1
        b1 = b1_biphasic_aqPhase85to1
        IW1_1_temp = IW1_1_temp_biphasic_aqPhase85to1

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_biphasic_aqPhase85to1
        LW2_1 = LW2_1_biphasic_aqPhase85to1

        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_biphasic_aqPhase85to1
        y1_step1%gain = y1_step1_gain_biphasic_aqPhase85to1
        y1_step1%xoffset = y1_step1_xoffset_biphasic_aqPhase85to1

        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1 ! samples

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_biphasic_aqPhase85to1_KGv4

!****************************************************************************************

    SUBROUTINE inverted_NNBAT_biphasic_aqPhase_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(20)             :: b1
        REAL(wrp), DIMENSION(60)             :: IW1_1_temp
        REAL(wrp), DIMENSION(20,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,20)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_biphasic_aqPhase
        x1_step1%gain = x1_step1_gain_biphasic_aqPhase
        x1_step1%ymin = x1_step1_ymin_biphasic_aqPhase

        ! Layer 1
        b1 = b1_biphasic_aqPhase
        IW1_1_temp = IW1_1_temp_biphasic_aqPhase

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_biphasic_aqPhase
        LW2_1 = LW2_1_biphasic_aqPhase

        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_biphasic_aqPhase
        y1_step1%gain = y1_step1_gain_biphasic_aqPhase
        y1_step1%xoffset = y1_step1_xoffset_biphasic_aqPhase

        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_biphasic_aqPhase_KGv4

!****************************************************************************************

    SUBROUTINE inverted_NNBAT_biphasic_orgPhase85to1_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(30)             :: b1
        REAL(wrp), DIMENSION(90)             :: IW1_1_temp
        REAL(wrp), DIMENSION(30,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,30)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_biphasic_orgPhase85to1
        x1_step1%gain = x1_step1_gain_biphasic_orgPhase85to1
        x1_step1%ymin = x1_step1_ymin_biphasic_orgPhase85to1

        ! Layer 1
        b1 = b1_biphasic_orgPhase85to1
        IW1_1_temp = IW1_1_temp_biphasic_orgPhase85to1

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
        IW1_1(i,:) = IW1_1_temp(m:n)
        m = m + SIZE(IW1_1, 2)
        n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_biphasic_orgPhase85to1
        LW2_1 = LW2_1_biphasic_orgPhase85to1

        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_biphasic_orgPhase85to1
        y1_step1%gain = y1_step1_gain_biphasic_orgPhase85to1
        y1_step1%xoffset = y1_step1_xoffset_biphasic_orgPhase85to1

        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1 ! samples

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_biphasic_orgPhase85to1_KGv4

!****************************************************************************************

    SUBROUTINE inverted_NNBAT_biphasic_orgPhase_KGv4(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(20)             :: b1
        REAL(wrp), DIMENSION(60)             :: IW1_1_temp
        REAL(wrp), DIMENSION(20,3)           :: IW1_1
        REAL(wrp), DIMENSION(1)              :: b2
        REAL(wrp), DIMENSION(1,20)           :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, repmaty1, &
                                                & repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(3), x1_step1%gain(3), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_biphasic_orgPhase
        x1_step1%gain = x1_step1_gain_biphasic_orgPhase
        x1_step1%ymin = x1_step1_ymin_biphasic_orgPhase

        ! Layer 1
        b1 = b1_biphasic_orgPhase
        IW1_1_temp = IW1_1_temp_biphasic_orgPhase

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_biphasic_orgPhase
        LW2_1 = LW2_1_biphasic_orgPhase

        ! Output 1
        ALLOCATE(y1_step1%xoffset(1), y1_step1%gain(1), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_biphasic_orgPhase
        y1_step1%gain = y1_step1_gain_biphasic_orgPhase
        y1_step1%xoffset = y1_step1_xoffset_biphasic_orgPhase

        ! ===== SIMULATION ========

        ! Dimensions
        Q = 1

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3(1)

    END SUBROUTINE inverted_NNBAT_biphasic_orgPhase_KGv4

!****************************************************************************************

    SUBROUTINE invert_NN_85to90_transfer_weights(a_w, transfert, weight_lower, weight_higher)

        IMPLICIT NONE
        
        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: a_w
        
        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: transfert, weight_lower, weight_higher

        ! Local variables
        REAL(wrp)                            :: upper_stop, lower_stop, normalize_x_value, &
                                                & zero_shift, scale_factor, x_span, min_val, &
                                                & max_val, scaled_sigmoid

        upper_stop = 0.90_wrp
        lower_stop = 0.85_wrp

        IF (a_w(1) < lower_stop) THEN
            transfert = [0.0_wrp]
            weight_lower = [1.0_wrp]
            weight_higher = [0.0_wrp]

        ELSEIF (a_w(1) > upper_stop) THEN
            transfert = [0.0_wrp]
            weight_lower = [0.0_wrp]
            weight_higher = [1.0_wrp]

        ELSE
            normalize_x_value = (a_w(1) - lower_stop)/(upper_stop - lower_stop)
            zero_shift = 0.5_wrp
            scale_factor = 10.0_wrp
            x_span = 0.5_wrp

            min_val = 0.00669285092428486_wrp ! 1/(1+exp(-scale_factor.*((zero_shift-x_span)-zero_shift))
            max_val = 0.986614298151430_wrp ! 1/(1+exp(-scale_factor.*((zero_shift+x_span)-zero_shift)))-min_val

            scaled_sigmoid = (1.0_wrp/(1.0_wrp + EXP( -scale_factor * (normalize_x_value - zero_shift ))) - min_val) / max_val

            transfert = [1.0_wrp]
            weight_lower = [1.0_wrp - scaled_sigmoid]
            weight_higher = [scaled_sigmoid]
        END IF

    END SUBROUTINE invert_NN_85to90_transfer_weights

!****************************************************************************************

    SUBROUTINE NN_VBSBAT_layers20_individProp_v2(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(20)             :: b1
        REAL(wrp), DIMENSION(900)            :: IW1_1_temp
        REAL(wrp), DIMENSION(20,45)          :: IW1_1
        REAL(wrp), DIMENSION(11)             :: b2
        REAL(wrp), DIMENSION(220)            :: LW2_1_temp
        REAL(wrp), DIMENSION(11,20)          :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(45), x1_step1%gain(45), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_layers20_individProp
        x1_step1%gain = x1_step1_gain_layers20_individProp
        x1_step1%ymin = x1_step1_ymin_layers20_individProp

        ! Layer 1
        b1 = b1_layers20_individProp
        IW1_1_temp = IW1_1_temp_layers20_individProp

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_layers20_individProp
        LW2_1_temp = LW2_1_temp_layers20_individProp

        m = 1
        n = SIZE(LW2_1, 2)
        DO i = 1, SIZE(LW2_1, 1)
            LW2_1(i,:) = LW2_1_temp(m:n)
            m = m + SIZE(LW2_1, 2)
            n = n + SIZE(LW2_1, 2)
        END DO

        ! Output 1
        ALLOCATE(y1_step1%xoffset(11), y1_step1%gain(11), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_layers20_individProp
        y1_step1%gain = y1_step1_gain_layers20_individProp

        y1_step1%xoffset = y1_step1_xoffset_layers20_individProp

        ! ===== SIMULATION ========

        ! Dimensions

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3

    END SUBROUTINE NN_VBSBAT_layers20_individProp_v2

!****************************************************************************************

    SUBROUTINE NN_VBSBAT_layers24_meanProp_v1(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(24)             :: b1
        REAL(wrp), DIMENSION(600)            :: IW1_1_temp
        REAL(wrp), DIMENSION(24,25)          :: IW1_1
        REAL(wrp), DIMENSION(11)             :: b2
        REAL(wrp), DIMENSION(264)            :: LW2_1_temp
        REAL(wrp), DIMENSION(11,24)          :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(25), x1_step1%gain(25), x1_step1%ymin(1))
        x1_step1%xoffset = x1_step1_xoffset_layers24_meanProp
        x1_step1%gain = x1_step1_gain_layers24_meanProp
        x1_step1%ymin = x1_step1_ymin_layers24_meanProp

        ! Layer 1
        b1 = b1_layers24_meanProp
        IW1_1_temp = IW1_1_temp_layers24_meanProp

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_layers24_meanProp
        LW2_1_temp = LW2_1_temp_layers24_meanProp

        m = 1
        n = SIZE(LW2_1, 2)
        DO i = 1, SIZE(LW2_1, 1)
            LW2_1(i,:) = LW2_1_temp(m:n)
            m = m + SIZE(LW2_1, 2)
            n = n + SIZE(LW2_1, 2)
        END DO

        ! Output 1
        ALLOCATE(y1_step1%xoffset(11), y1_step1%gain(11), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_layers24_meanProp
        y1_step1%gain = y1_step1_gain_layers24_meanProp
        y1_step1%xoffset = y1_step1_xoffset_layers24_meanProp


        ! ===== SIMULATION ========

        ! Dimensions

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
            repmatx2(m:n) = x1_step1%gain
            m = m + SIZE(x1_step1%gain, 1)
            n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
            repmatx3(m:n) = x1_step1%ymin
            m = m + SIZE(x1_step1%ymin, 1)
            n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3

    END SUBROUTINE NN_VBSBAT_layers24_meanProp_v1

!****************************************************************************************

    SUBROUTINE NN_VBSBAT_layers28_individProp_v1(x1, y1)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: x1

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: y1

        ! Local variables
        TYPE(step_input)                     :: x1_step1, y1_step1
        REAL(wrp), DIMENSION(28)             :: b1
        REAL(wrp), DIMENSION(1260)           :: IW1_1_temp
        REAL(wrp), DIMENSION(28,45)          :: IW1_1
        REAL(wrp), DIMENSION(11)             :: b2
        REAL(wrp), DIMENSION(308)            :: LW2_1_temp
        REAL(wrp), DIMENSION(11,28)          :: LW2_1
        REAL(wrp), DIMENSION(:), ALLOCATABLE :: xp1_1, xp1_2, xp1, repmatx1, repmatx2, repmatx3, y1_1, y1_2, &
                                                & repmaty1, repmaty2, repmaty3, a1, a2
        INTEGER(wip)                         :: i, n, m, Q, p

        ! Input 1
        ALLOCATE(x1_step1%xoffset(45), x1_step1%gain(45), x1_step1%ymin(1))

        x1_step1%xoffset = x1_step1_xoffset_layers28_individProp
        x1_step1%gain = x1_step1_gain_layers28_individProp
        x1_step1%ymin = x1_step1_ymin_layers28_individProp

        ! Layer 1
        b1 = b1_layers28_individProp
        IW1_1_temp = IW1_1_temp_layers28_individProp

        m = 1
        n = SIZE(IW1_1, 2)
        DO i = 1, SIZE(IW1_1, 1)
            IW1_1(i,:) = IW1_1_temp(m:n)
            m = m + SIZE(IW1_1, 2)
            n = n + SIZE(IW1_1, 2)
        END DO

        ! Layer 2
        b2 = b2_layers28_individProp
        LW2_1_temp = LW2_1_temp_layers28_individProp

        m = 1
        n = SIZE(LW2_1, 2)
        DO i = 1, SIZE(LW2_1, 1)
            LW2_1(i,:) = LW2_1_temp(m:n)
            m = m + SIZE(LW2_1, 2)
            n = n + SIZE(LW2_1, 2)
        END DO

        ! Output 1
        ALLOCATE(y1_step1%xoffset(11), y1_step1%gain(11), y1_step1%ymin(1))
        y1_step1%ymin = y1_step1_ymin_layers28_individProp
        y1_step1%gain = y1_step1_gain_layers28_individProp
        y1_step1%xoffset =  y1_step1_xoffset_layers28_individProp


        ! ===== SIMULATION ========

        ! Dimensions

        ! Input 1
        p = SIZE(x1_step1%gain, 1)/ SIZE(x1, 1)
        ! Step 1
        ALLOCATE(repmatx1(SIZE(x1_step1%xoffset)*p))
        repmatx1 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%xoffset, 1)
        DO i = 1, p
            repmatx1(m:n) = x1_step1%xoffset
            m = m + SIZE(x1_step1%xoffset, 1)
            n = n + SIZE(x1_step1%xoffset, 1)
        END DO
        ALLOCATE(xp1_1(SIZE(x1, 1)))
        xp1_1 = x1 - repmatx1

        ! Step 2
        ALLOCATE(repmatx2(SIZE(x1_step1%gain)*p))
        repmatx2 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%gain, 1)
        DO i = 1, p
        repmatx2(m:n) = x1_step1%gain
        m = m + SIZE(x1_step1%gain, 1)
        n = n + SIZE(x1_step1%gain, 1)
        END DO
        ALLOCATE(xp1_2(SIZE(x1, 1)))
        xp1_2 = xp1_1 * repmatx2

        ! Step 3
        ALLOCATE(repmatx3(SIZE(x1_step1%ymin)*p))
        repmatx3 = 0.0_wrp
        m = 1
        n = SIZE(x1_step1%ymin, 1)
        DO i = 1, p
        repmatx3(m:n) = x1_step1%ymin
        m = m + SIZE(x1_step1%ymin, 1)
        n = n + SIZE(x1_step1%ymin, 1)
        END DO
        ALLOCATE(xp1(SIZE(x1, 1)))
        xp1 = xp1_2 + repmatx3(1)

        ! Layer 1
        ALLOCATE(a1(SIZE(b1, 1)))
        a1 = 2.0_wrp/(1.0_wrp + EXP(-2.0_wrp * (b1 + MATMUL(IW1_1,xp1)))) - 1.0_wrp

        ! Layer 2
        ALLOCATE(a2(1))
        a2 = b2 + MATMUL(LW2_1,a1)

        ! Output 1
        p = SIZE(y1_step1%gain, 1)/ SIZE(a2, 1)

        ! Step 1
        ALLOCATE(repmaty1(SIZE(y1_step1%ymin)*p))
        repmaty1 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%ymin, 1)
        DO i = 1, p
            repmaty1(m:n) = y1_step1%ymin
            m = m + SIZE(y1_step1%ymin, 1)
            n = n + SIZE(y1_step1%ymin, 1)
        END DO
        ALLOCATE(y1_1(SIZE(a2, 1)))
        y1_1 = a2 - repmaty1(1)

        ! Step 2
        ALLOCATE(repmaty2(SIZE(y1_step1%gain)*p))
        repmaty2 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%gain, 1)
        DO i = 1, p
            repmaty2(m:n) = y1_step1%gain
            m = m + SIZE(y1_step1%gain, 1)
            n = n + SIZE(y1_step1%gain, 1)
        END DO
        ALLOCATE(y1_2(SIZE(a2, 1)))
        y1_2 = y1_1/repmaty2

        ! Step 3
        ALLOCATE(repmaty3(SIZE(y1_step1%xoffset)*p))
        repmaty3 = 0.0_wrp
        m = 1
        n = SIZE(y1_step1%xoffset, 1)
        DO i = 1, p
            repmaty3(m:n) = y1_step1%xoffset
            m = m + SIZE(y1_step1%xoffset, 1)
            n = n + SIZE(y1_step1%xoffset, 1)
        END DO

        y1 = y1_2 + repmaty3

    END SUBROUTINE NN_VBSBAT_layers28_individProp_v1

!****************************************************************************************

END MODULE NN_MOD
