!****************************************************************************************
!*   MODULE OBJECTIVE_FUNCTION_MOD                                                      *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that defines the objective function, which is based on the error for        *
!*   convergence of the volatility basis set (VBS) (i.e., based on the error of the     *
!*   fraction of each organic species that is partitioned to the condensed phase (Ej)). *
!*   The objective function is the nonlinear function which is used by the equation     *
!*   solver Powell hybrid method (hybrd1) to find a 0.                                  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

MODULE OBJECTIVE_FUNCTION_MOD

    USE PRECISION_MOD ! Module that defines the precision/length of REAL/INTERGER/CHAR data type
	USE TOOLS_MOD ! Module that defines molar mass of water and other useful subroutines/functions
    USE OBJECTIVE_FUNCTION_VARIABLES_MOD ! Module that defines the variables of the objective function that
                                         ! are used/updated at each iteration

    IMPLICIT NONE

CONTAINS

    !****************************************************************************************
    !*   SUBROUTINE error_out_powell_Csat												    *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that needs to be provided to the solver of nonlinear equations hybrd1   *
    !*   (Powell hybrid method). This subroutine calculates the functions.                  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************
        
    SUBROUTINE error_out_powell_Csat(n, guess_Ej_objfun, error_out_objfun, iflag)

        IMPLICIT NONE
        
        ! Input
        INTEGER(4), INTENT(IN)               :: n                ! Number of organic species in the system
                                                                 ! The number of functions and variables in the system of equations
        ! Output
        INTEGER(4), INTENT(OUT)              :: iflag            ! Output flag related to function evaluation
        REAL(8), DIMENSION(n), INTENT(OUT)   :: error_out_objfun ! Error based on the fraction of species that partitions to the condensed phase (Ej)
                                                                 ! The function value evaluated at the ouput guess_Ej_objfun
        ! Input/Output
        REAL(8), DIMENSION(n), INTENT(INOUT) :: guess_Ej_objfun  ! The fraction of org. species j that partitions to the condensed phase (function variable)
                                                                 ! Input: initial estimate of the solution vector
                                                                 ! Output: final estimate of the solution vector
                                  

        ! Local variables
        REAL(8), DIMENSION(n)              :: mass_fraction_water_alpha, mass_fraction_water_beta, &
                                                & activity_coefficient_alpha, activity_coefficient_beta, &
                                                & q_alpha, q_beta, mass_fraction_water_alpha_denominator, &
                                                & mass_fraction_water_beta_denominator, alpha_good_denominators, &
                                                & not_alpha_good_denominators, beta_good_denominators, &
                                                & not_beta_good_denominators, Coa_j, Coa_j_alpha, Coa_j_beta, &
                                                & Cstar_j_objfun_via_alpha, Cstar_j_objfun_via_beta, Cstar_j_objfun_used, &
                                                & Ej_new, Coa_j_new, Coa_j_alpha_new, Caq_j_alpha_new, Coa_j_beta_new, Caq_j_beta_new
        REAL(8), DIMENSION(1)              :: Coa_guess_viaEj, Caq_alpha, Coaaq_alpha, Caq_beta, Coaaq_beta, Coaaq_via_Ej_guess, &
                                                & massweighted_molar_weight_alpha, massweighted_molar_weight_beta, Coa_new_viaEj, &
                                                & Caq_alpha_new, Coa_alpha_new, Caq_beta_new, Coa_beta_new, q_alpha_water_new
    
        iflag = 0

        ! Make sure the bounds of guess_Ej_objfun are under control (between 0.0 and 1.0)
        !CALL EjBoundsCheck_Hard(guess_Ej_objfun)

		guess_Ej_objfun = EjBoundsCheck_Soft(guess_Ej_objfun, 0.0_wrp, 1.0_wrp)

        ! Extract terms for phases alpha and beta from previous iteration
        activity_coefficient_alpha = activity_coefficient_AB_objfun(1,:) ! Activity coefficient of each organic species in phase alpha
        activity_coefficient_beta = activity_coefficient_AB_objfun(2,:) ! Activity coefficient of each organic species in phase beta

        q_alpha = q_alpha_molefrac_phase_split_org_objfun ! Fraction of organic species j in phase alpha (relative to total of j in all liquid phases)
        q_beta = 1.0_wrp - q_alpha ! Fraction of organic species j in phase beta (relative to total of j in all liquid phases)

        mass_fraction_water_alpha = mass_fraction_water_AB_objfun(1,:) ! Mass fraction of water in phase alpha (associated with each organic species)
        mass_fraction_water_beta = mass_fraction_water_AB_objfun(2,:) ! Mass fraction of water in phase beta (associated with each organic species)
        mass_fraction_water_alpha_denominator = (1.0_wrp - mass_fraction_water_alpha) ! Mass fraction of each organic species j in phase alpha
        mass_fraction_water_beta_denominator = (1.0_wrp - mass_fraction_water_beta) ! Mass fraction of each organic species j in phase beta

        ! Check that denominator is not negative. If denominator is negative (mass fraction organic species j > 1.0_wrp), it becomes equal to 1.0_wrp
        !alpha_good_denominators = 0.0_wrp
        !not_alpha_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_alpha_denominator > 0.0_wrp) alpha_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_alpha_denominator > 0.0_wrp) not_alpha_good_denominators = 0.0_wrp
        !beta_good_denominators = 0.0_wrp
        !not_beta_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_beta_denominator > 0.0_wrp) beta_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_beta_denominator > 0.0_wrp) not_beta_good_denominators = 0.0_wrp
        WHERE (mass_fraction_water_alpha_denominator < tinynumber) mass_fraction_water_alpha_denominator = 1.0_wrp
        WHERE (mass_fraction_water_beta_denominator < tinynumber) mass_fraction_water_beta_denominator = 1.0_wrp
		
        !! Apply correction to the mass fraction of organic species j (if needed)
        !mass_fraction_water_alpha_denominator = ABS(mass_fraction_water_alpha_denominator) * alpha_good_denominators + not_alpha_good_denominators
        !mass_fraction_water_beta_denominator = ABS(mass_fraction_water_beta_denominator) * beta_good_denominators + not_beta_good_denominators

        ! New Coa_objfun_j calculation
        Coa_j = guess_Ej_objfun * C_OM_ugPm3_objfun ! Mass concentration of each organic species j that is in the condensed phase
                                                    ! C_particle_j = Ej * C_particle+gas_j = Ej * C_total_j
        Coa_guess_viaEj(1) = SUM(Coa_j, 1) ! Total organic mass concentration of all organics species that are in the condensed phase
                                           ! Cparticle_allspecies = SUM(Cparticle_j)

        ! Calculate the water mass concentration in phase alpha associated with each organic species j
        Coa_j_alpha = Coa_j * q_alpha ! Mass concentration in phase alpha of each organic species j
                                      ! Cparticle-alpha_j = Cparticle_j * q-alpha_j

        ! mass_fraction_water_alpha_j / mass_fraction_water_alpha_denominator_j = mass_water_alpha_j/mass_org_alpha_j
        Caq_alpha(1) = SUM(Coa_j_alpha * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator, 1)

        ! Calculate the water mass concentration in phase beta associated with each organic species j
        Coa_j_beta = Coa_j * q_beta ! organic mass in phase beta of species j
                                    ! Cparticle-beta_j = Cparticle_j * q-beta_j

        ! mass_fraction_water_beta_j / mass_fraction_water_beta_denominator_j = mass_water_beta_j/mass_org_beta_j
        Caq_beta(1) = SUM(Coa_j_beta * mass_fraction_water_beta / mass_fraction_water_beta_denominator, 1)

        ! Calculate total mass concentration (phase + beta; water + organic)
        Coaaq_alpha(1) = SUM(Coa_j_alpha, 1) + Caq_alpha(1) ! Total mass concentraiton in phase alpha: water mass + organic mass
        Coaaq_beta(1) = SUM(Coa_j_beta, DIM = 1) + Caq_beta(1) ! Total mass concentraiton in phase beta: water mass + organic mass
        Coaaq_via_Ej_guess(1) = Coaaq_beta(1) + Coaaq_alpha(1) ! Total mass concentration (phase + beta; water + organic)
                                                               ! Coaaq_via_Ej_guess is calculated from input Ej

        ! Calculate the effective saturation concentration of each species C*_j via phase alpha
        massweighted_molar_weight_alpha(1) = SUM((Coa_j_alpha/(molecular_weight_objfun)), 1) + Caq_alpha(1)/M_water

        IF (massweighted_molar_weight_alpha(1) > 0.0_wrp) THEN
            ! See equation 7 of paper
            Cstar_j_objfun_via_alpha = Cstar_dry_objfun * activity_coefficient_alpha * q_alpha !* &
				!& Coaaq_via_Ej_guess(1)/(molecular_weight_objfun * massweighted_molar_weight_alpha(1))
        ELSE
            Cstar_j_objfun_via_alpha = 0.0_wrp
        END IF

        ! Calculate the effective saturation concentration of each species C*_j via phase beta
        massweighted_molar_weight_beta(1) = SUM((Coa_j_beta/(molecular_weight_objfun)), 1) + Caq_beta(1)/M_water

        IF (massweighted_molar_weight_beta(1) > 0.0_wrp) THEN
            ! See equation 7 of paper
            Cstar_j_objfun_via_beta = Cstar_dry_objfun * activity_coefficient_beta * q_beta !* &
				!& Coaaq_via_Ej_guess(1)/(molecular_weight_objfun * massweighted_molar_weight_beta(1))
        ELSE
            Cstar_j_objfun_via_beta = 0.0_wrp
        END IF

        ! Select C* to use (weighted average by q_alpha)
        Cstar_j_objfun_used = Cstar_j_objfun_via_alpha * q_alpha + Cstar_j_objfun_via_beta * q_beta

        ! new Ej calculation (See equation 8 of paper)
        Ej_new = 1.0_wrp/(1.0_wrp + Cstar_j_objfun_used/Coaaq_via_Ej_guess(1))

        ! Calculate new mass concentration values based on new Ej value.
        Coa_j_new = Ej_new * C_OM_ugPm3_objfun ! Cparticle_j = Ej * Cparticle+gas_j
        Coa_new_viaEj(1) = SUM(Coa_j_new, DIM = 1) ! Total organic mass of all species that is in the condensed phase
                                                   ! Cparticle_allspecies = SUM(Cparticle_j)
        ! New alpha mass concentrations
        Coa_j_alpha_new = Coa_j * q_alpha ! Mass concentration in phase alpha of each organic species j
        Caq_j_alpha_new = (Coa_j_alpha_new * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator) ! Now for water
        Coa_alpha_new(1) = SUM(Coa_j_alpha_new, 1) ! Total organic mass concentration in phase alpha
        Caq_alpha_new(1) = SUM(Caq_j_alpha_new, 1) ! Total water mass concentration in phase alpha

        ! New beta mass concentrations
        Coa_j_beta_new = Coa_j_new * q_beta ! Mass concentration in phase beta of each organic species j
        Caq_j_beta_new = (Coa_j_beta_new * mass_fraction_water_beta / mass_fraction_water_beta_denominator) ! Now for water
        Coa_beta_new(1) = SUM(Coa_j_beta_new, 1) ! Total organic mass concentration in phase beta
        Caq_beta_new(1) = SUM(Caq_j_beta_new, 1) ! Total water mass concentration in phase beta
        
        ! Calculate the new fraction of water for each organic species j in liquid phase alpha
        IF (Caq_alpha_new(1) + Caq_beta_new(1) < tinynumber) THEN
            q_alpha_water_new(1) = 0.0_wrp
        ELSE
            q_alpha_water_new(1) = Caq_alpha_new(1)/(Caq_alpha_new(1) + Caq_beta_new(1))
		END IF
        
        ! Calculate relative error between old Ej and new Ej in terms of logarithms.
		! Make sure numerator and denominator are not equal to 0.0_wrp
        WHERE (guess_Ej_objfun(1:n) < 0.0_wrp + tinynumber .OR. Ej_new(1:n) < 0.0_wrp + tinynumber) ! == 0.0_wrp
            error_out_objfun(1:n) = LOG((guess_Ej_objfun(1:n)+tinynumber)/(Ej_new(1:n)+tinynumber))
        ELSEWHERE
            error_out_objfun(1:n) = LOG(guess_Ej_objfun(1:n)/Ej_new(1:n)) ! System of equations
		ENDWHERE
		
        !error_out_objfun(1:n) = (guess_Ej_objfun(1:n) - Ej_new(1:n))**2 
    
	END SUBROUTINE error_out_powell_Csat
	
	!****************************************************************************************
    !*   SUBROUTINE error_out_powell												        *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Subroutine that needs to be provided to the solver of nonlinear equations hybrd1   *
    !*   (Powell hybrid method). This subroutine calculates the functions.                  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************
        
    SUBROUTINE error_out_powell(n, guess_Ej_objfun, error_out_objfun, iflag)

        IMPLICIT NONE
        
        ! Input
        INTEGER(4), INTENT(IN)               :: n                ! Number of organic species in the system
                                                                 ! The number of functions and variables in the system of equations
        ! Output
        INTEGER(4), INTENT(OUT)              :: iflag            ! Output flag related to function evaluation
        REAL(8), DIMENSION(n), INTENT(OUT)   :: error_out_objfun ! Error based on the fraction of species that partitions to the condensed phase (Ej)
                                                                 ! The function value evaluated at the ouput guess_Ej_objfun
        ! Input/Output
        REAL(8), DIMENSION(n), INTENT(INOUT) :: guess_Ej_objfun  ! The fraction of org. species j that partitions to the condensed phase (function variable)
                                                                 ! Input: initial estimate of the solution vector
                                                                 ! Output: final estimate of the solution vector
                                  

        ! Local variables
        REAL(8), DIMENSION(n)              :: mass_fraction_water_alpha, mass_fraction_water_beta, &
                                                & activity_coefficient_alpha, activity_coefficient_beta, &
                                                & q_alpha, q_beta, mass_fraction_water_alpha_denominator, &
                                                & mass_fraction_water_beta_denominator, alpha_good_denominators, &
                                                & not_alpha_good_denominators, beta_good_denominators, &
                                                & not_beta_good_denominators, Coa_j, Coa_j_alpha, Coa_j_beta, &
                                                & Cstar_j_objfun_via_alpha, Cstar_j_objfun_via_beta, Cstar_j_objfun_used, &
                                                & Ej_new, Coa_j_new, Coa_j_alpha_new, Caq_j_alpha_new, Coa_j_beta_new, Caq_j_beta_new
        REAL(8), DIMENSION(1)              :: Coa_guess_viaEj, Caq_alpha, Coaaq_alpha, Caq_beta, Coaaq_beta, Coaaq_via_Ej_guess, &
                                                & massweighted_molar_weight_alpha, massweighted_molar_weight_beta, Coa_new_viaEj, &
                                                & Caq_alpha_new, Coa_alpha_new, Caq_beta_new, Coa_beta_new, q_alpha_water_new
    
        iflag = 0

        ! Make sure the bounds of guess_Ej_objfun are under control (between 0.0 and 1.0)
        CALL EjBoundsCheck_Hard(guess_Ej_objfun)
		!guess_Ej_objfun = EjBoundsCheck_Soft(guess_Ej_objfun, 0.0_wrp, 1.0_wrp)

        ! Extract terms for phases alpha and beta from previous iteration
        activity_coefficient_alpha = activity_coefficient_AB_objfun(1,:) ! Activity coefficient of each organic species in phase alpha
        activity_coefficient_beta = activity_coefficient_AB_objfun(2,:) ! Activity coefficient of each organic species in phase beta

        q_alpha = q_alpha_molefrac_phase_split_org_objfun ! Fraction of organic species j in phase alpha (relative to total of j in all liquid phases)
        q_beta = 1.0_wrp - q_alpha ! Fraction of organic species j in phase beta (relative to total of j in all liquid phases)

        mass_fraction_water_alpha = mass_fraction_water_AB_objfun(1,:) ! Mass fraction of water in phase alpha (associated with each organic species)
        mass_fraction_water_beta = mass_fraction_water_AB_objfun(2,:) ! Mass fraction of water in phase beta (associated with each organic species)
        mass_fraction_water_alpha_denominator = (1.0_wrp - mass_fraction_water_alpha) ! Mass fraction of each organic species j in phase alpha
        mass_fraction_water_beta_denominator = (1.0_wrp - mass_fraction_water_beta) ! Mass fraction of each organic species j in phase beta

        ! Check that denominator is not negative. If denominator is negative (mass fraction organic species j > 1.0_wrp), it becomes equal to 1.0_wrp
        !alpha_good_denominators = 0.0_wrp
        !not_alpha_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_alpha_denominator > 0.0_wrp) alpha_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_alpha_denominator > 0.0_wrp) not_alpha_good_denominators = 0.0_wrp
        !beta_good_denominators = 0.0_wrp
        !not_beta_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_beta_denominator > 0.0_wrp) beta_good_denominators = 1.0_wrp
        !WHERE (mass_fraction_water_beta_denominator > 0.0_wrp) not_beta_good_denominators = 0.0_wrp
        WHERE (mass_fraction_water_alpha_denominator < tinynumber) mass_fraction_water_alpha_denominator = 1.0_wrp
        WHERE (mass_fraction_water_beta_denominator < tinynumber) mass_fraction_water_beta_denominator = 1.0_wrp
		
        !! Apply correction to the mass fraction of organic species j (if needed)
        !mass_fraction_water_alpha_denominator = ABS(mass_fraction_water_alpha_denominator) * alpha_good_denominators + not_alpha_good_denominators
        !mass_fraction_water_beta_denominator = ABS(mass_fraction_water_beta_denominator) * beta_good_denominators + not_beta_good_denominators

        ! New Coa_objfun_j calculation
        Coa_j = guess_Ej_objfun * C_OM_ugPm3_objfun ! Mass concentration of each organic species j that is in the condensed phase
                                                    ! C_particle_j = Ej * C_particle+gas_j = Ej * C_total_j
        Coa_guess_viaEj(1) = SUM(Coa_j, 1) ! Total organic mass concentration of all organics species that are in the condensed phase
                                           ! Cparticle_allspecies = SUM(Cparticle_j)

        ! Calculate the water mass concentration in phase alpha associated with each organic species j
        Coa_j_alpha = Coa_j * q_alpha ! Mass concentration in phase alpha of each organic species j
                                      ! Cparticle-alpha_j = Cparticle_j * q-alpha_j

        ! mass_fraction_water_alpha_j / mass_fraction_water_alpha_denominator_j = mass_water_alpha_j/mass_org_alpha_j
        Caq_alpha(1) = SUM(Coa_j_alpha * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator, 1)

        ! Calculate the water mass concentration in phase beta associated with each organic species j
        Coa_j_beta = Coa_j * q_beta ! organic mass in phase beta of species j
                                    ! Cparticle-beta_j = Cparticle_j * q-beta_j

        ! mass_fraction_water_beta_j / mass_fraction_water_beta_denominator_j = mass_water_beta_j/mass_org_beta_j
        Caq_beta(1) = SUM(Coa_j_beta * mass_fraction_water_beta / mass_fraction_water_beta_denominator, 1)

        ! Calculate total mass concentration (phase + beta; water + organic)
        Coaaq_alpha(1) = SUM(Coa_j_alpha, 1) + Caq_alpha(1) ! Total mass concentraiton in phase alpha: water mass + organic mass
        Coaaq_beta(1) = SUM(Coa_j_beta, DIM = 1) + Caq_beta(1) ! Total mass concentraiton in phase beta: water mass + organic mass
        Coaaq_via_Ej_guess(1) = Coaaq_beta(1) + Coaaq_alpha(1) ! Total mass concentration (phase + beta; water + organic)
                                                               ! Coaaq_via_Ej_guess is calculated from input Ej

        ! Calculate the effective saturation concentration of each species C*_j via phase alpha
        massweighted_molar_weight_alpha(1) = SUM((Coa_j_alpha/(molecular_weight_objfun)), 1) + Caq_alpha(1)/M_water

        IF (massweighted_molar_weight_alpha(1) > 0.0_wrp) THEN
            ! See equation 7 of paper
            Cstar_j_objfun_via_alpha = Cstar_dry_objfun * activity_coefficient_alpha * q_alpha * &
				& Coaaq_via_Ej_guess(1)/(molecular_weight_objfun * massweighted_molar_weight_alpha(1))
        ELSE
            Cstar_j_objfun_via_alpha = 0.0_wrp
        END IF

        ! Calculate the effective saturation concentration of each species C*_j via phase beta
        massweighted_molar_weight_beta(1) = SUM((Coa_j_beta/(molecular_weight_objfun)), 1) + Caq_beta(1)/M_water

        IF (massweighted_molar_weight_beta(1) > 0.0_wrp) THEN
            ! See equation 7 of paper
            Cstar_j_objfun_via_beta = Cstar_dry_objfun * activity_coefficient_beta * q_beta * &
				& Coaaq_via_Ej_guess(1)/(molecular_weight_objfun * massweighted_molar_weight_beta(1))
        ELSE
            Cstar_j_objfun_via_beta = 0.0_wrp
        END IF

        ! Select C* to use (weighted average by q_alpha)
        Cstar_j_objfun_used = Cstar_j_objfun_via_alpha * q_alpha + Cstar_j_objfun_via_beta * q_beta

        ! new Ej calculation (See equation 8 of paper)
        Ej_new = 1.0_wrp/(1.0_wrp + Cstar_j_objfun_used/Coaaq_via_Ej_guess(1))

        ! Calculate new mass concentration values based on new Ej value.
        Coa_j_new = Ej_new * C_OM_ugPm3_objfun ! Cparticle_j = Ej * Cparticle+gas_j
        Coa_new_viaEj(1) = SUM(Coa_j_new, DIM = 1) ! Total organic mass of all species that is in the condensed phase
                                                   ! Cparticle_allspecies = SUM(Cparticle_j)
        ! New alpha mass concentrations
        Coa_j_alpha_new = Coa_j * q_alpha ! Mass concentration in phase alpha of each organic species j
        Caq_j_alpha_new = (Coa_j_alpha_new * mass_fraction_water_alpha / mass_fraction_water_alpha_denominator) ! Now for water
        Coa_alpha_new(1) = SUM(Coa_j_alpha_new, 1) ! Total organic mass concentration in phase alpha
        Caq_alpha_new(1) = SUM(Caq_j_alpha_new, 1) ! Total water mass concentration in phase alpha

        ! New beta mass concentrations
        Coa_j_beta_new = Coa_j_new * q_beta ! Mass concentration in phase beta of each organic species j
        Caq_j_beta_new = (Coa_j_beta_new * mass_fraction_water_beta / mass_fraction_water_beta_denominator) ! Now for water
        Coa_beta_new(1) = SUM(Coa_j_beta_new, 1) ! Total organic mass concentration in phase beta
        Caq_beta_new(1) = SUM(Caq_j_beta_new, 1) ! Total water mass concentration in phase beta
        
        ! Calculate the new fraction of water for each organic species j in liquid phase alpha
        IF (Caq_alpha_new(1) + Caq_beta_new(1) < tinynumber) THEN
            q_alpha_water_new(1) = 0.0_wrp
        ELSE
            q_alpha_water_new(1) = Caq_alpha_new(1)/(Caq_alpha_new(1) + Caq_beta_new(1))
		END IF
        
        ! Calculate relative error between old Ej and new Ej in terms of logarithms.
		! Make sure numerator and denominator are not equal to 0.0_wrp
        WHERE (guess_Ej_objfun(1:n) < 0.0_wrp + tinynumber .OR. Ej_new(1:n) < 0.0_wrp + tinynumber) ! == 0.0_wrp
            error_out_objfun(1:n) = LOG((guess_Ej_objfun(1:n)+tinynumber)/(Ej_new(1:n)+tinynumber))
        ELSEWHERE
            error_out_objfun(1:n) = LOG(guess_Ej_objfun(1:n)/Ej_new(1:n)) ! System of equations
		ENDWHERE
		
        !error_out_objfun(1:n) = (guess_Ej_objfun(1:n) - Ej_new(1:n))**2 
    
    END SUBROUTINE error_out_powell

    !****************************************************************************************
    !*   SUBROUTINE EjBoundsCheck                                                           *
    !*                                                                                      *
    !*   :: Purpose ::                                                                      *
    !*   Makes sure that the input to Powell's method is a number between 0.0 and 1.0.      *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
    !*   McGill University, Montreal, Quebec (2021),                                        *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created: 15-04-2021                                                             *
    !*   -> latest changes:                                                                 *
    !****************************************************************************************

    PURE ELEMENTAL SUBROUTINE EjBoundsCheck_Hard(Ej_set)

        IMPLICIT NONE

        ! Input/Output
        REAL(8), INTENT(INOUT) :: Ej_set                ! Fraction of each organic species that is
                                                          ! partitioned to the condensed phase
        ! Local variables
        REAL(8), PARAMETER     :: scaleval = 1.0E-4_wrp ! Small value used to adjust lower/upper bounds of Ej_set

        ! First use smart adjustments on values with lower/upper bounds if necessary:
        IF (Ej_set < 0.0_wrp) THEN
            Ej_set = ABS(Ej_set)*scaleval
            IF (Ej_set > 1.0_wrp) THEN
                Ej_set = MIN(Ej_set*tinynumber, scaleval)
            ENDIF
        ELSE IF (Ej_set > 1.0_wrp) THEN
            Ej_set = 1.0_wrp -(Ej_set-1.0_wrp)*scaleval
            IF (Ej_set < 0.0_wrp) THEN
                Ej_set = MAX(ABS(Ej_set)*tinynumber, 1.0_wrp-scaleval)
            ENDIF
        ENDIF
            
	END SUBROUTINE EjBoundsCheck_Hard

    !********************************************************************************************
    !*   :: Purpose ::                                                                          *
    !*  Elemental function to constrain an input value x to fall within provided lower and      *
    !*  upper limits. If x exceeds one of the bounds, the value of x will be adjusted in a      *
    !*  responsive manner, such that subsequent exceeding, yet different values of x would lead *
    !*  to a different output (not simply setting hard bounds). Hard limits are only enforced   *
    !*  in cases in which the "softness" of the applied bounds is exceeded by extreme x values. *
    !*                                                                                          *
    !*   :: Author & Copyright ::                                                               *
    !*   Andi Zuend,                                                                            *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                              *
    !*                                                                                          *
    !*   -> created:        2022-11-04                                                          *
    !*   -> latest changes: 2022-11-04                                                          *
    !*                                                                                          *
    !********************************************************************************************
    PURE ELEMENTAL FUNCTION EjBoundsCheck_Soft(x, l_lim, u_lim) result(y)
    
    IMPLICIT NONE
    !interface arguments:
    REAL(wrp),INTENT(IN) :: x            !input; considered a natural log-scale value    
    REAL(wrp),INTENT(IN) :: l_lim        !imposed hard lower limit
    REAL(wrp),INTENT(IN) :: u_lim        !imposed hard upper limit
    REAL(wrp)            :: y            !output; x value, but potentially truncated to meet set bound constraints
    !local variables:
    REAL(wrp),parameter  :: scaler = 1.0E-3_wrp 
    REAL(wrp)            :: local_low_lim, local_upper_lim, margin, trunc
    !.............................
    
    !(1) determine margin of x values in which soft bounds apply:
	! l_lim	= epsilon(0.0_wrp)
	! u_lim = 1.0_wrp - l_lim
    margin = 1.0E-3_wrp*ABS(u_lim - l_lim)
    local_low_lim = l_lim + margin
    local_upper_lim = u_lim - margin
    
    !(2) apply soft bounds:
    IF (x > local_upper_lim) THEN
        trunc = scaler*LOG(ABS(x))
        y = local_upper_lim + trunc     !apply a scaled limit
        IF (y > u_lim) THEN             !very rare case; apply hard bounds
            y = u_lim    
        ENDIF
    ELSE IF (x < local_low_lim) THEN
        trunc = scaler*LOG(ABS(x))
        y = local_low_lim - trunc
        IF (y < l_lim) THEN
            y = l_lim    
        ENDIF 
    ELSE !no bound constraint necessary
        y = x
    END IF                         
    
    END FUNCTION EjBoundsCheck_Soft

END MODULE OBJECTIVE_FUNCTION_MOD
