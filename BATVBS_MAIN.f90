!****************************************************************************************
!*   FILE BATVBS.f90                                                                    *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Contains the Fortran subroutines that are used as interface to create Python       *
!*   functions through F2PY (https://numpy.org/doc/stable/f2py/).                       *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano Damha                                                               *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 01-10-2021                                                      *
!****************************************************************************************

!****************************************************************************************
!*   SUBROUTINE BAT                                                                     *
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
!*   -> latest changes: 01-10-2021                                                      *
!****************************************************************************************

SUBROUTINE BAT(M_org, O2C_org, H2C_org, N2C_org, Group_org, x_org, y_water, y_org)

	USE BATVBS_MOD      ! Module that contains the main subroutines to perform BAT and VBS calculations
	IMPLICIT NONE

    INTEGER, PARAMETER      :: ip = KIND(1.0E0)             ! Precision of INTEGERS                
    INTEGER, PARAMETER      :: rp = KIND(1.0D0)             ! Precision of REALS
	! inputs
    INTEGER(ip), INTENT(IN) :: Group_org                    ! Type of oxygen-bearing functional group that represents the organic species
    REAL(rp), INTENT(IN)    :: O2C_org, H2C_org, N2C_org, & ! Elemental ratios (O/C, H/C, N/C) of the organic species in the binary mixture  
                                & M_org, x_org                ! Molecular weight (g/mol) and mole fraction of the organic species in the binary mixture
	! outputs
    REAL(rp), INTENT(OUT)   :: y_water, y_org               ! Activity cofficient of water and organic species in the binary mixture
	
	! variables used to map precision of REALS
	REAL(wrp)               :: O2C_org_wrp, H2C_org_wrp, N2C_org_wrp, x_org_wrp, M_org_wrp, y_org_wrp, y_water_wrp
        
	! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
	O2C_org_wrp = REAL(O2C_org, kind = wrp)
	H2C_org_wrp = REAL(H2C_org, kind = wrp)
	N2C_org_wrp = REAL(N2C_org, kind = wrp)
	x_org_wrp = REAL(x_org, kind = wrp)
	M_org_wrp = REAL(M_org, kind = wrp)
	y_org_wrp = REAL(y_org, kind = wrp)
	y_water_wrp = REAL(y_water, kind = wrp)

	! calculates the activity cofficient of water and organic species in the binary mixture (org+water):  y_water_wrp, y_org_wrp
	CALL BAT_SIMULATION(M_org_wrp, O2C_org_wrp, H2C_org_wrp, N2C_org_wrp, Group_org, x_org_wrp, y_water_wrp, y_org_wrp)
        
	! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
	y_org = REAL(y_org_wrp, kind = rp)
	y_water = REAL(y_water_wrp, kind = rp)

END SUBROUTINE BAT
 
!****************************************************************************************
!*   SUBROUTINE BAT_INVERSION                                                           *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Subroutine that calculates the (mole-fraction-based) activity coefficient of each  *
!*   organic species j in liquid phases alpha (water-rich) and beta (organic-rich),     *
!*   the mass fraction of each organic species j in liquid phases alpha (water-rich)    *
!*   and beta (organic-rich), and the total (organic + water) mass concentration in     *
!*   liquid phases alpha and beta (ug/m3) at a given water activity (RH) using the      *
!*   BAT (INVERSION) model.                                                             *
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
!*   -> latest changes: 01-10-2021                                                      *
!****************************************************************************************

SUBROUTINE BAT_INVERSION(N_org, M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, Group_org_j, C_PM_org_j, a_water, &
        & y_org_j_alpha, y_org_j_beta, w_org_j_alpha, w_org_j_beta, C_PM_water_org_alpha, C_PM_water_org_beta)

	USE BATVBS_MOD	! Module that contains the subroutines to perform BAT and VBS calculations
    IMPLICIT NONE

    INTEGER, PARAMETER  						:: ip = KIND(1.0E0)                				! Precision of INTEGERS                
	INTEGER, PARAMETER  						:: rp = KIND(1.0D0)                 			! Precision of REALS
    ! Inputs
    INTEGER(ip), INTENT(IN)                     :: N_org                                     	! Number of organic species in the system
    INTEGER(ip), DIMENSION(N_org), INTENT(IN)   :: Group_org_j                                  ! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(rp), DIMENSION(N_org), INTENT(IN)      :: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, &  ! Molecular weight (g/mol) and elemental ratios (O/C, H/C, N/C) of each organic species j in the system
													& C_PM_org_j                                ! Total "dry" mass concentration of each organic species j in the particle (ug/m3)
    REAL(rp), INTENT(IN)                        :: a_water                                      ! Target relative humidity (gas-phase water activity) that the mole fraction of water (and organics) in the mixture must match (scale of 0 to 1)

    ! Outputs
    REAL(rp), INTENT(OUT)                       :: C_PM_water_org_alpha, C_PM_water_org_beta    ! Total (organic + water) mass concentration in liquid phases alpha and beta (ug/m3)
    REAL(rp), DIMENSION(N_org), INTENT(OUT)     :: y_org_j_alpha, y_org_j_beta, &               ! Activity coefficient of organic species j in liquid phases alpha and beta
													& w_org_j_alpha, w_org_j_beta               ! Mass fraction of organic species j in liquid phases alpha and beta
    
			
    
    ! variables used to map precision of REALS
    REAL(wrp)                                   :: a_water_wrp, C_PM_water_org_alpha_wrp, C_PM_water_org_beta_wrp
    REAL(wrp), DIMENSION(N_org)                	:: O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, &
                                                    & M_org_j_wrp, C_PM_org_j_wrp, &
                                                    & y_org_j_alpha_wrp, y_org_j_beta_wrp, &
                                                    & w_org_j_alpha_wrp, w_org_j_beta_wrp
    
	! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
	O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
	H2C_org_j_wrp = REAL(H2C_org_j, kind = wrp)
	N2C_org_j_wrp = REAL(N2C_org_j, kind = wrp)
	M_org_j_wrp = REAL(M_org_j, kind = wrp)
	C_PM_org_j_wrp = REAL(C_PM_org_j, kind = wrp)
	y_org_j_alpha_wrp = REAL(y_org_j_alpha, kind = wrp)
	y_org_j_beta_wrp = REAL(y_org_j_beta, kind = wrp)   
	w_org_j_alpha_wrp = REAL(w_org_j_alpha, kind = wrp)
	w_org_j_beta_wrp = REAL(w_org_j_beta, kind = wrp)
	a_water_wrp = REAL(a_water, kind = wrp)
	C_PM_water_org_alpha_wrp = REAL(C_PM_water_org_alpha, kind = wrp)
	C_PM_water_org_beta_wrp = REAL(C_PM_water_org_beta, kind = wrp)
			
	! calculates activity coefficient and mass fraction of each organic species j in liquid phases alpha (water-rich) and beta (water-poor): y_org_j_alpha, y_org_j_beta, w_org_j_alpha, w_org_j_beta
	! calculates the total (all organics + water) mass concentration in liquid phases alpha and beta: C_PM_water_org_alpha, C_PM_water_org_beta
	CALL BAT_INVERSION_SIMULATION(N_org, M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, Group_org_j, &
		& C_PM_org_j_wrp, a_water_wrp, y_org_j_alpha_wrp, y_org_j_beta_wrp, w_org_j_alpha_wrp, w_org_j_beta_wrp, &
		& C_PM_water_org_alpha_wrp, C_PM_water_org_beta_wrp)
		
	! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
	y_org_j_alpha = REAL(y_org_j_alpha_wrp, kind = rp)
	y_org_j_beta = REAL(y_org_j_beta_wrp, kind = rp)   
	w_org_j_alpha = REAL(w_org_j_alpha_wrp, kind = rp)
	w_org_j_beta = REAL(w_org_j_beta_wrp, kind = rp)
	C_PM_water_org_alpha = REAL(C_PM_water_org_alpha_wrp, kind = rp)
	C_PM_water_org_beta = REAL(C_PM_water_org_beta_wrp, kind = rp)
		
END SUBROUTINE BAT_INVERSION
    
!****************************************************************************************
!*   SUBROUTINE BAT_INVERSION_VBS                                                       *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Subroutine that calculates the equilibrium aerosol mass concentration (ug/m3) of   *
!*   water and each organic species j in phases alpha (water-rich) and beta (water-poor)*
!*   at a given water activity (RH) using the coupled BAT(INVERSION)+VBS model.         *
!*   (Gorkowski, K., Preston, T. C., and Zuend, A.: Relative-humidity-dependent organic *
!*   aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem.     *
!*   Phys., 19, 13383–13407, https://doi.org/10.5194/acp-19-13383-2019, 2019.).         *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 01-10-2021                                                      *
!****************************************************************************************

SUBROUTINE BAT_INVERSION_VBS(N_a_water, N_org, M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, C_sat_org_j, C_gas_particle_org_j, &
    & Group_org_j, a_water, C_PM_org_j_alpha, C_PM_org_j_beta, C_PM_water_alpha, C_PM_water_beta)

    USE BATVBS_MOD ! Module that contains the subroutines to perform VBS+BAT calculation
    IMPLICIT NONE

	INTEGER, PARAMETER  						:: ip = KIND(1.0E0)                				! Precision of INTEGERS                
	INTEGER, PARAMETER  						:: rp = KIND(1.0D0)                 			! Precision of REALS
    ! inputs
	INTEGER(ip), INTENT(IN)						:: N_a_water									! Number of water activities 
    INTEGER(ip), INTENT(IN)                     :: N_org                                    	! Number of organic species j in the system
    INTEGER(ip), DIMENSION(N_org), INTENT(IN)   :: Group_org_j									! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(rp), DIMENSION(N_org), INTENT(IN)      :: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, &	! Molar mass (g/mol) and elemental ratios (O/C, H/C) of each organic species j in the system
													& C_gas_particle_org_j, &                	! Mass concentration of each organic species j in gas phase + particle phase (ug/m3)
                                                    & C_sat_org_j                            	! Pure-component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
    REAL(rp), DIMENSION(N_a_water), INTENT(IN)	:: a_water                                  	! Target relative humidity (gas-phase water activity) that the mole fraction of water (and organics) in the mixture must match

    ! output
    REAL(rp), DIMENSION(N_a_water), INTENT(OUT)				:: C_PM_water_alpha, C_PM_water_beta      		! Equilibrium water mass concentration in liquid phases alpha and beta that results from VBS calculations (ug/m3)
    REAL(rp), DIMENSION(N_a_water, N_org), INTENT(OUT)		:: C_PM_org_j_alpha, C_PM_org_j_beta      		! Equilibrium mass concentration of each organic species j in liquid phases alpha and beta that results from VBS calculations (ug/m3)
        		 
    ! variables used to map precision of REALS
    REAL(wrp), DIMENSION(N_a_water)				:: C_PM_water_alpha_wrp, C_PM_water_beta_wrp, &
													& a_water_wrp
    REAL(wrp), DIMENSION(N_org)                 :: O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, &
													& M_org_j_wrp, C_gas_particle_org_j_wrp, &
													& C_sat_org_j_wrp
	
	REAL(wrp), DIMENSION(N_a_water, N_org)		:: C_PM_org_j_alpha_wrp, C_PM_org_j_beta_wrp
	
	! Write to file the actual inputs we get from python
	
    ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
    O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
    H2C_org_j_wrp = REAL(H2C_org_j, kind = wrp)
    N2C_org_j_wrp = REAL(N2C_org_j, kind = wrp)
    M_org_j_wrp = REAL(M_org_j, kind = wrp)
    C_gas_particle_org_j_wrp = REAL(C_gas_particle_org_j, kind = wrp)
    C_sat_org_j_wrp = REAL(C_sat_org_j, kind = wrp)
    a_water_wrp = REAL(a_water, kind = wrp)
    C_PM_org_j_alpha_wrp = REAL(C_PM_org_j_alpha, kind = wrp)
    C_PM_org_j_beta_wrp = REAL(C_PM_org_j_beta, kind = wrp)
    C_PM_water_alpha_wrp = REAL(C_PM_water_alpha, kind = wrp)
    C_PM_water_beta_wrp = REAL(C_PM_water_beta, kind = wrp)
	
    ! calculates the total equilibrium OA mass concentration of organics and water in phases alpha and beta	
  !  CALL BAT_INVERSION_VBS_SIMULATION(N_org, M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, &
		!& C_sat_org_j_wrp, C_gas_particle_org_j_wrp, Group_org_j, a_water_wrp, C_PM_org_j_alpha_wrp, &
		!& C_PM_org_j_beta_wrp, C_PM_water_alpha_wrp, C_PM_water_beta_wrp)
    CALL BAT_INVERSION_VBS_SIMULATION(N_a_water, N_org, M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, &
		& C_sat_org_j_wrp, C_gas_particle_org_j_wrp, Group_org_j, a_water_wrp, C_PM_org_j_alpha_wrp, &
		& C_PM_org_j_beta_wrp, C_PM_water_alpha_wrp, C_PM_water_beta_wrp)
        
    ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
    C_PM_org_j_alpha = REAL(C_PM_org_j_alpha_wrp, kind = rp)
    C_PM_org_j_beta = REAL(C_PM_org_j_beta_wrp, kind = rp)
    C_PM_water_alpha = REAL(C_PM_water_alpha_wrp, kind = rp)
    C_PM_water_beta = REAL(C_PM_water_beta_wrp, kind = rp)
	
END SUBROUTINE BAT_INVERSION_VBS
	
!****************************************************************************************
!*   SUBROUTINE CSAT_APPROX					                                            *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Estimates the pure component saturation mass concentration of an organic species   *
!*	 j at the temperature of interest, C_sat_org_j (ug/m3), at dry conditions (RH = 0)	*
!*   when C_star_org_j (ug/m3) and C_PM_org_j (ug/m3) are known.						*
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 19-02-2022                                                             *
!*   -> latest changes: 23-02-2022                                                      *
!****************************************************************************************

SUBROUTINE CSAT_APPROX(N_org, M_org_j, O2C_org_j, C_star_org_j, C_gas_particle_org_j, Group_org_j, C_sat_org_j)

    USE BATVBS_MOD ! Module that contains the subroutines to perform VBS+BAT calculation
    IMPLICIT NONE

	!INTEGER, PARAMETER  						:: ip =	KIND(1.0E0)                				! Precision of INTEGERS                
	!INTEGER, PARAMETER  						:: rp = KIND(1.0D0)                 			! Precision of REALS
    ! inputs
    INTEGER(wip), INTENT(IN)                     :: N_org                                    	! Number of organic species j in the system
    INTEGER(wip), DIMENSION(N_org), INTENT(IN)   :: Group_org_j									! Type of oxygen-bearing functional group of each organic species j in the system
    REAL(wrp), DIMENSION(N_org), INTENT(IN)      :: M_org_j, O2C_org_j, & 	                    ! Molar mass (g/mol) and O/C elemental ratio of each organic species j in the system
													& C_gas_particle_org_j, &                	! Mass concentration of each organic species j in gas phase + particle phase (ug/m3)
                                                    & C_star_org_j                            	! Effective concentration C* of each organic species j at the temperature of interest (ug/m3)


    ! output
    REAL(wrp), DIMENSION(N_org), INTENT(OUT)		:: C_sat_org_j                            	    ! Pure component saturation mass concentration of each organic species j at the temperature of interest (ug/m3)
    ! variables used to map precision of REALS
    REAL(wrp), DIMENSION(N_org)                 :: O2C_org_j_wrp, M_org_j_wrp, C_gas_particle_org_j_wrp, &
													& C_star_org_j_wrp, C_sat_org_j_wrp
												
    ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
 !   O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
 !   M_org_j_wrp = REAL(M_org_j, kind = wrp)
 !   C_gas_particle_org_j_wrp = REAL(C_gas_particle_org_j, kind = wrp)
	!C_star_org_j_wrp = REAL(C_star_org_j, kind = wrp)
 !   C_sat_org_j_wrp = REAL(C_sat_org_j, kind = wrp)
	
    ! calculates the total equilibrium OA mass concentration of organics and water in phases alpha and beta
	!CALL CSAT_ESTIMATION(N_org, M_org_j_wrp, O2C_org_j_wrp, C_star_org_j_wrp, C_gas_particle_org_j_wrp, &
	!	& Group_org_j, C_sat_org_j_wrp)
		CALL CSAT_ESTIMATION(N_org, M_org_j, O2C_org_j, C_star_org_j, C_gas_particle_org_j, &
		& Group_org_j, C_sat_org_j)
        
    ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
    !C_sat_org_j = REAL(C_sat_org_j_wrp, kind = rp)
	
END SUBROUTINE CSAT_APPROX

!****************************************************************************************
!*   SUBROUTINE DENSITY_ORG_ESTIMATE                                               *
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

SUBROUTINE DENSITY_ORG_ESTIMATE(N_org, M_org_j, O2C_org_j, H2C_org_j, N2C_org_j, density_org_j)

    USE BATVBS_MOD ! Module that contains the subroutines to perform VBS+BAT calculation
	IMPLICIT NONE
	
	INTEGER, PARAMETER  						:: ip = KIND(1.0E0)                				! Precision of INTEGERS                
	INTEGER, PARAMETER  						:: rp = KIND(1.0D0)                 			! Precision of REALS
	! inputs
	INTEGER(ip), INTENT(IN) 					:: N_org										! Number of organic species j in the system
	REAL(rp), DIMENSION(N_org), INTENT(IN)  	:: M_org_j, O2C_org_j, H2C_org_j, N2C_org_j 	! Elemental ratios (O/C, H/C, N/C)
																								! and molecular weight of each org species j

	! outputs
	REAL(rp), DIMENSION(N_org), INTENT(OUT)		:: density_org_j       							! Density estimate of each org species j

    ! variables used to map precision of REALS
    REAL(wrp), DIMENSION(N_org)                 	:: O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, &
													& M_org_j_wrp, density_org_j_wrp
													
    ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
    O2C_org_j_wrp = REAL(O2C_org_j, kind = wrp)
    H2C_org_j_wrp = REAL(H2C_org_j, kind = wrp)
    N2C_org_j_wrp = REAL(N2C_org_j, kind = wrp)
    M_org_j_wrp = REAL(M_org_j, kind = wrp)
	density_org_j_wrp = REAL(density_org_j, kind = wrp)
	
    ! estimates the density of each organic compound j
    CALL Org_density_Estimate_KGv1(M_org_j_wrp, O2C_org_j_wrp, H2C_org_j_wrp, N2C_org_j_wrp, density_org_j_wrp)
	
    ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
    density_org_j = REAL(density_org_j_wrp, kind = rp)

END SUBROUTINE DENSITY_ORG_ESTIMATE