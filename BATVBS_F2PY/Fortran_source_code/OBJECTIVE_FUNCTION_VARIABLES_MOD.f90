!****************************************************************************************
!*   MODULE OBJECTIVE_FUNCTION_VARIABLES_MOD                                            *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that declares the variables of the objective function (error for            *
!*   convergence of volatility basis set (VBS), based on error of the fraction of each  *
!*   organic species that is partitioned to the condensed phase (Ej)).                  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

MODULE OBJECTIVE_FUNCTION_VARIABLES_MOD

    USE PRECISION_MOD ! Module that defines the precision/length of REAL/INTERGER/CHAR data type

    IMPLICIT NONE

    REAL(wrp), DIMENSION(:), ALLOCATABLE   :: C_OM_ugPm3_objfun, &                         ! Total mass concentration of organics (particle+gas)
                                              & Cstar_dry_objfun, &                        ! Saturation mass concensation of organics
                                              & q_alpha_molefrac_phase_split_org_objfun, & ! Fraction of organic species j in phase alpha
                                              & molecular_weight_objfun                    ! Molecular weight of each org. j

    REAL(wrp), DIMENSION(:,:), ALLOCATABLE :: activity_coefficient_AB_objfun, &            ! Activity coefficient of org. j in phases alpha and beta
                                              mass_fraction_water_AB_objfun                ! Mass fraction of water in phases alpha and beta for each org. j

    REAL(wrp), DIMENSION(:), ALLOCATABLE   :: partition_coefficients_objfun, &             ! Fraction of org. j partitioned to the condensed phase (all liquids)
                                              & Cstar_j_objfun, &                          ! Effective saturation concentration of each org. j
                                              & Coa_AB_objfun, &                           ! Total organic mass concentration in phases alpha and beta
                                              & Caq_AB_objfun, &                           ! Total water mass concentration in phases alpha and beta
                                              & Coa_objfun, &                              ! Total organic mass concentration (phase alpha + phase beta)
                                              & q_alpha_water_objfun                       ! Fraction of water in phase alpha for each org. j

    REAL(wrp), DIMENSION(:,:), ALLOCATABLE :: Coa_j_AB_objfun, &                           ! Mass concentration of each org. j in phases alpha and beta
                                              & Caq_j_AB_objfun                            ! Mass concentration of water in phases alpha and beta for each org.j
 
    REAL(wrp), DIMENSION(:), ALLOCATABLE   :: error_out_objfun                             ! Error based on the fraction of species that partitions to the                                                                                         ! condensed phase (Ej). The function value evaluated at the ouput                                                                                       ! guess_Ej_objfun

	REAL(wrp), DIMENSION(:), ALLOCATABLE   :: partition_coefficients_memory_effect		   ! Partitioning coefficient variable used to test memory effect
	INTEGER(wip)						   :: flag_memory_effect						   ! VBS solver status flag
END MODULE OBJECTIVE_FUNCTION_VARIABLES_MOD


    
