!****************************************************************************************
!*   MODULE TOOLS_MOD                                                                   *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that contains subroutines to perform different recurrent calculations       *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

MODULE TOOLS_MOD

    USE PRECISION_MOD ! Module that defines the precision/length of REAL/INTERGER/CHAR data type

    IMPLICIT NONE
    
    REAL(wrp), PARAMETER, PUBLIC  :: M_water = 18.015_wrp ! Molecular weight of water

CONTAINS

    !****************************************************************************************
    !*   SUBROUTINE interp1                                                                 *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that performs 1-D data interpolation                                    *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE interp1(xin, yin, xout, yout)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: xin, yin, xout ! xin vector contains the sample points,
                                                              ! yin vector contains the corresponding values yin(xin),
                                                              ! xout vector contains the coordinates of the query points.
        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: yout ! yout vector contains the corresponding values yout(xout)

        ! Local variables
        REAL(wrp), DIMENSION(SIZE(xin, 1))   :: xin_new, yin_new
        REAL(wrp) :: TEMP_X, TEMP_Y
        INTEGER(wip) :: inputIndex, L, K

        IF (xout(1) > MAXVAL(xin, 1) .OR. xout(1) < MINVAL(xin, 1)) THEN
            yout(1) = -99999.0_wrp
        ELSE IF (xout(1) <= MAXVAL(xin, 1) .AND. xout(1) >= MINVAL(xin, 1)) THEN
            xin_new = xin
            yin_new = yin
            ! Sort in increasing order
            DO K = 1, SIZE(xin,1) - 1
                DO L = K + 1, SIZE(xin,1)
                    IF (xin_new(K).GT.xin_new(L)) THEN
                        TEMP_X = xin_new(K)
                        TEMP_Y = yin_new(K)
                        xin_new(K) = xin_new(L)
                        yin_new(K) = yin_new(L)
                        xin_new(L) = TEMP_X
                        yin_new(L) = TEMP_Y
                    END IF
                END DO
            END DO

            DO inputIndex = 1, (SIZE(xin,1)-1)
                IF (xout(1) >= xin_new(inputIndex) .AND. xout(1) < xin_new(inputIndex+1)) THEN
                    yout(1) = yin_new(inputIndex) + (xout(1) - xin_new(inputIndex))*(yin_new(inputIndex+1) - &
                        yin_new(inputIndex))/(xin_new(inputIndex+1) - xin_new(inputIndex))
                END IF
            END DO

        END IF

    END SUBROUTINE interp1

    !****************************************************************************************
    !*   SUBROUTINE molar_based_means                                                       *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Calculates the molar and mass weighted average properties of a system from the     *
    !*   individual component properties.                                                   *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************
    
    SUBROUTINE molar_based_means(C_ug_m3, M_gPmol, O2C, H2C, N2C, &
        & avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C)

        IMPLICIT NONE
        
        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: O2C, H2C, N2C, & ! Elemental ratios (O/C, H/C, N/C) of org.
                                                & C_ug_m3, M_gPmol ! Mass concentration and molecular weight of org.
        ! Outputs
        REAL(wrp), INTENT(OUT)               :: avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C
                                                ! Average molecular weight; average O/C, H/C and H/C elemental ratios
        ! Local variables
        REAL(wrp), DIMENSION(SIZE(M_gPmol ,1))     :: H2Cest, NC, NO, NH, NN, molar_C_ug_m3, mass_fraction_inPM, N2Cest
        REAL(wrp) :: MC, MO, MH, MN, molar_NC
        INTEGER(wip) :: sizeM
        
        sizeM = SIZE(M_gPmol)

        ! Molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]
        MC = 12.0100_wrp
        MO = 16.0000_wrp
        MH = 1.00800_wrp
        MN = 14.0067_wrp
        
        H2Cest = H2C
        N2Cest = N2C

        ! 1) Set the N2C value if not provided from input (if a negative value is given at input)
        WHERE (N2C < 0.0_wrp) N2Cest = 0.0_wrp
        ! 2) Estimate the H2C value if not provided from input
        WHERE (H2C < 0.0_wrp) H2Cest = 2.0_wrp - O2C ! estimate H2C assuming an aliphatic compound with H2C = 2 in absence of oxygen-bearing
                                                        ! functional groups, then correct for oxygen content assuming a -1 slope (Van Krevelen Diagram of typical SOA).

        ! 2) Compute the approx. number of atoms per organic molecule
        NC = M_gPmol/(MC + H2Cest*MH + O2C*MO + N2Cest*MN) ! Number of carbon
        NO = O2C*NC ! Number of oxygen
        NH = H2Cest*NC ! Number of hydrogen
        NN = N2Cest*NC ! Number of nitrogen
        
        ! 3) Compute molar concentration
        molar_C_ug_m3 = C_ug_m3/M_gPmol ! mol/m3
        
        ! 4) New ratios based on molar concentrations
        molar_NC = SUM(molar_C_ug_m3*NC)
        molar_avg_O2C = SUM(molar_C_ug_m3*NO)/molar_NC
        molar_avg_H2C = SUM(molar_C_ug_m3*NH)/molar_NC
        molar_avg_N2C = SUM(molar_C_ug_m3*NN)/molar_NC
        
        ! 5) Compute average molar weight
        mass_fraction_inPM = C_ug_m3/SUM(C_ug_m3)
        avg_M_gPmol = SUM(mass_fraction_inPM * M_gPmol)
        mass_weighted_avg_O2C = SUM(mass_fraction_inPM * O2C)

    END SUBROUTINE molar_based_means
    
    !****************************************************************************************
    !*   SUBROUTINE Org_density_Estimate_KGv1                                               *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that computes the density of pure liquid organic                        *
    !*   compounds at 298.15 K using a simple structure-based method                        *
    !*   (G. S. Girolami J. Chem. Educ. 71, 962-964, 1994).                                 *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************
 
    SUBROUTINE Org_density_Estimate_KGv1(M, O2C, H2C, N2C, densityEst)

        IMPLICIT NONE
        
        ! Inputs
        REAL(wrp), DIMENSION(:), INTENT(IN)  :: M, O2C, H2C, N2C ! elemental ratios (O/C, H/C, N/C)
                                                                 ! and molecular weight of organic species

        ! Outputs
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: densityEst       ! density estimate for each organic species

        ! Local variables
        REAL(wrp), DIMENSION(SIZE(M ,1))     :: H2Cest, NC, rho1, N2Cest
        REAL(wrp) :: MC, MO, MH, MN
        INTEGER(wip) :: i, sizeM

        sizeM = SIZE(M,1)

        !  DensityEst(M, O2C, H2C) ,  H:C input is optional (when known, otherwise enter a negative value)
        ! the molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]
        MC = 12.0100_wrp
        MO = 16.0000_wrp
        MH = 1.00800_wrp
        MN = 14.0067_wrp
        
        N2Cest = N2C
        H2Cest = H2C
        
        ! 1) Set the N2C value if not provided from input (if a negative value is given at input)
        WHERE (N2C < 0.0_wrp) N2Cest = 0.0_wrp
        
        ! 2) Estimate the H2C value if not provided from input (if a negative value is given at input)
        WHERE (H2C < 0.0_wrp) H2Cest = 2.0_wrp - O2C ! estimate H2C assuming an aliphatic compound with H2C = 2 in absence of oxygen-bearing.
                                                         ! functional groups, then correct for oxygen content assuming a-1 slope (Van Krevelen Diagram of typical SOA).

        ! 3) Compute the approx. number of carbon atoms per organic molecule
        NC = M/(MC + H2Cest*MH + O2C*MO + N2Cest*MN)

        ! 4) compute density estimate based on method by Girolami (1994)
        !  here no correction is applied for rings and aromatic compounds (due to limited info at input)
        rho1 = M/(5.0_wrp*NC*(2.0_wrp + H2Cest + O2C*2.0_wrp + N2Cest*2.0_wrp))

        DO i = 1, SIZE(NC, 1)
            densityEst(i) = rho1(i)*(1.000_wrp + MINVAL([NC(i) * O2C(i) * 0.100_wrp + NC(i) * N2Cest(i) * &
                0.100_wrp, 0.300_wrp], 1))
        ! density in [g/cm^3]; here it is scaled assuming that most of the oxygen atoms are able
        ! to make H-bonds (donor or acceptor).
        END DO

    END SUBROUTINE Org_density_Estimate_KGv1

!****************************************************************************************
!*   SUBROUTINE replace_data_A_to_B_KGv1                                                *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Replaces an element of an array by another element specified at input.             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

    SUBROUTINE replace_data_A_to_B_KGv1(y, z, w, v)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN) :: y ! input array
        REAL(wrp), INTENT(IN)               :: z ! element to be replaced
        REAL(wrp), INTENT(IN)               :: w ! replacement value

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT) :: v ! resulting array

        ! Local variables
        INTEGER(wip) :: i

        v = y

        DO i = 1, SIZE(y, 1)
            IF (y(i) == z) THEN
                v(i) = w
            END IF
        END DO

    END SUBROUTINE replace_data_A_to_B_KGv1

    !****************************************************************************************
    !*   SUBROUTINE convert_chemical_structure_to_OH_eqv_v3                                 *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Translate the properties of functionalized molecules to a hypothetical             *
    !*   OH-equivalent molecule of modified O/C and molecular weight using shift fit data.  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE convert_chemical_structure_to_OH_eqv_v3(O2C, molarmass_ratio, shift_method, O2C_eqv, molarmass_ratio_eqv)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT(IN)       :: O2C, molarmass_ratio ! Elemental O/C ratio and molecular weight of each org. species
        CHARACTER(clen), DIMENSION(:), INTENT(IN) :: shift_method ! functional group of each organic species

        ! Output
        REAL(wrp), DIMENSION(:), INTENT(OUT)      :: O2C_eqv, molarmass_ratio_eqv ! Translated molecule properties (to OH-equivalent):
                                                                                  ! Hypothetical OH-equivalent O/C and molecular weight ratios

        ! Local variables
        REAL(wrp), DIMENSION(4)                   :: fit_shift
        REAL(wrp)                                 :: MW_old, MW_new
        INTEGER(wip)                              :: sp_i, calc_shift

        DO sp_i = 1, SIZE(O2C, 1)
            calc_shift = 1_wip

            IF  ((shift_method(sp_i) == 'hydroxyl') .OR. (shift_method(sp_i) == 'carboxyl')) THEN
                calc_shift = 0
            ELSE IF ((shift_method(sp_i) == 'hydroperoxideSOA') .OR. (shift_method(sp_i) == 'SOA chemicals')) THEN
                fit_shift = [0.000149021455554774_wrp, 0.00473627706154738_wrp, 0.869057801919811_wrp, 0.564783492434737_wrp]
            ELSE IF (shift_method(sp_i) == 'PEG') THEN
                fit_shift = [0.00544768078879267_wrp, 3.86433605822883_wrp, -0.267168022244528_wrp, 0.255486696379870_wrp]
            ELSE IF (shift_method(sp_i) == 'ketone') THEN
                fit_shift = [0.00453425820008037_wrp, 0.000648450309348991_wrp, 0.138144408029760_wrp, 0.352454367906330_wrp]
            ELSE IF (shift_method(sp_i) == 'hydroperoxide') THEN
                fit_shift = [8.17160348517250E-06_wrp, 4.53182593743281E-07_wrp, 0.966089559236154_wrp, 0.459433193460024_wrp]
            ELSE IF (shift_method(sp_i) == 'ether') THEN
                fit_shift = [2.44336043644347E-05_wrp, 0.000158316167479487_wrp, 0.284974095922167_wrp, 0.229338647030993_wrp]
            ELSE IF (shift_method(sp_i) == 'ester') THEN
                fit_shift = [-1.29324633442629_wrp, 0.00108128314380665_wrp, 1.24051435479678_wrp, 0.405354156019572_wrp]
            ELSE
                PRINT *, 'no O2C and MW system selected'
                calc_shift = 0_wip
            END IF

            IF (calc_shift == 1) THEN
                O2C_eqv(sp_i) = O2C(sp_i)/(1.0_wrp + fit_shift(3)*EXP(-(O2C(sp_i))*fit_shift(1)))
                MW_old = (18.016_wrp/molarmass_ratio(sp_i))
                MW_new = MW_old/(1.0_wrp + fit_shift(4)*EXP(-(MW_old)*fit_shift(2)))
                molarmass_ratio_eqv(sp_i) = 18.016_wrp/MW_new
            ELSE
                O2C_eqv(sp_i) = O2C(sp_i)
                molarmass_ratio_eqv(sp_i) = molarmass_ratio(sp_i)
            END IF
        END DO

    END SUBROUTINE convert_chemical_structure_to_OH_eqv_v3

    !****************************************************************************************
    !*   SUBROUTINE replace_data_A_to_B_KGv1                                                *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Calculates exp(x) where x is within the interval [-690.7755_wrp,+690.7755_wrp]     *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    SUBROUTINE exp_wlimiter(y, v)

        IMPLICIT NONE

        ! Input
        REAL(wrp), DIMENSION(:), INTENT (IN)  :: y ! argument of exp(x)

        ! Output
        REAL(wrp), DIMENSION(:), INTENT (OUT) :: v ! value of exp(x)

        ! Local variables
        REAL(wrp)                             :: z, w
        INTEGER(wip)                          :: i

        ! Define an interval
        z = -690.7755_wrp
        w = 690.7755_wrp

        DO i = 1, SIZE(y, 1)
            IF (y(i) < z) THEN
                v(i) = z
            ELSE IF (y(i) > w) THEN
                v(i) = w
            ELSE
                v(i) = y(i)
            END IF
        END DO
        v = exp(v)

    END SUBROUTINE exp_wlimiter

    !------------------------------------------------------------------------------------------------------
    !** subroutine to determine existence and array indices of local minimum and maximum values of a given array **
    !** note: array must contain more than 3 points for this procedure to make sense. Also, only one local min & one max will be determined. **
    pure subroutine DetermineLocalMinMax(curvedat, minmaxfound, idx_min, idx_max)
    
        implicit none
        !interface arguments:
        real(wrp),dimension(:),intent(in) :: curvedat             !input: array containing curve data to be evalulated;
        logical(wlp),intent(out)          :: minmaxfound          !output: if .true. then curve contains both a local min and max; 
        integer(wip),intent(out)          :: idx_min, idx_max     !output: determined index of array element at approximate local minimum / maximum; set to 0 if not applicable;
        !local variables:
        integer(wip) :: i
        logical(wlp) :: adiff_test1, adiff_test2, locmin_found, locmax_found
        !.......................................

        !Approach: determine the indices of the limits of an interval containing a local minimum by means of
        !comparing the activity value at an index location to those to the "left" and "right" of the point.
        locmin_found = .false.
        locmax_found = .false.
        minmaxfound = .false.
        idx_min = 0
        idx_max = 0
        adiff_test1 = curvedat(2) < curvedat(1)
        do i = 2, size(curvedat)-1
            adiff_test2 = curvedat(i) < curvedat(i+1)
            if (adiff_test1) then
                if (adiff_test2) then                       !if .true. --> local minimum detected
                    locmin_found = .true.
                    idx_min = i
                endif
            else                                            !activity at index 'i' is greater than at previous array index
                if (.not. adiff_test2) then                 !for local maximum detection
                    locmax_found = .true.
                    idx_max = i
                endif
            endif
            if (locmin_found .and. locmax_found) then
                minmaxfound = .true.
                exit                                        !done searching the data array and found both local min and max; no reason to continue loop;
            endif
            adiff_test1 = .not. adiff_test2                 !transfer difference test result for re-use in next iteration;
        enddo

    end subroutine DetermineLocalMinMax
    !------------------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------------------
    !** function for computing the linear interpolation value of an x-value given the y-values (activities here) and a y-target value. **
    pure function LinearInterpolation(xval, yval, target_yval) result(intp_xval)
    
        implicit none
        real(wrp),dimension(2),intent(in) :: xval, yval
        real(wrp),intent(in) :: target_yval
        real(wrp) :: intp_xval
        !.........................
        intp_xval = ( xval(1)*(yval(1) - target_yval) + xval(2)*(target_yval - yval(2)) ) &
                & / (yval(1) - yval(2))
    
    end function LinearInterpolation
    !------------------------------------------------------------------------------------------------------

    subroutine linspace(from_value, to_value, array)
        real(wrp), intent(in) :: from_value, to_value
        real(wrp), intent(out) :: array(:)
        real(wrp) :: range
        integer(wip) :: n, i
        n = size(array)
        range = to_value - from_value

        if (n == 0) return

        if (n == 1) then
            array(1) = from_value
            return
        end if

        do i=1, n
            array(i) = from_value + range * (i - 1) / (n - 1)
        end do
    end subroutine
    
    function wtime ( )

    !*****************************************************************************80
    !
    !! WTIME returns a reading of the wall clock time.
    !
    !  Discussion:
    !
    !    To get the elapsed wall clock time, call WTIME before and after a given
    !    operation, and subtract the first reading from the second.
    !
    !    This function is meant to suggest the similar routines:
    !
    !      "omp_get_wtime ( )" in OpenMP,
    !      "MPI_Wtime ( )" in MPI,
    !      and "tic" and "toc" in MATLAB.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    27 April 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, real ( kind = wrp ) WTIME, the wall clock reading, in seconds.
    !
      implicit none

      integer (wip) clock_max
      integer (wip) clock_rate
      integer (wip) clock_reading
      real (wrp) wtime

      call system_clock ( clock_reading, clock_rate, clock_max )

      wtime = real ( clock_reading, kind = wrp ) &
            / real ( clock_rate, kind = wrp )

      return
    end

END MODULE TOOLS_MOD

