!****************************************************************************************
!*   PROGRAM EXAMPLE_BAT                                                                *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Fortran program that shows how to use the outputs of the BAT_SIMULATION subroutine.*
!*                                                                                      *
!*   SUBROUTINE BAT_SIMULATION                                                          *
!*   Input  : a binary mixture composed of 1 organic species and its associated         *
!*           water amount.                                                              * 
!*                                                                                      *
!*   Output : activity coefficients of water and organic species assuming both are      *
!*   present in the same single liquid phase (unlike more sophisticated BAT             *
!*   calculations with phase separation consideration).                                 *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 2021-04-15                                                             *
!*   -> latest changes: 2021-07-10                                                      *
!****************************************************************************************

PROGRAM EXAMPLE_BAT

    IMPLICIT NONE 
    
    !On variable KIND: For compatibility with typical choices in other Fortran programs and f2py,  
    !we use INTEGER(ip) and REAL(rp), where ip = KIND(1.0E0) and rp = KIND(1.0D0), as precision KINDs in this main program (serving as interface). 
    !Within the BAT procedures, a set working real precision (wrp) and similar for integer is used, 
    !which are defined in module PRECISION_MOD (and accessible via BATVBS_MOD).
    
    ! Inputs
    INTEGER, PARAMETER  :: ip = KIND(1.0E0)                 ! Precision of INTEGERS                
    INTEGER, PARAMETER  :: rp = KIND(1.0D0)                 ! Precision of REALS
    INTEGER(ip)         :: Group_org                        ! Type of oxygen-bearing functional group that represents the organic species (1 species) in the binary mixture (org+water)
    REAL(rp)            :: O2C_org, H2C_org, N2C_org, &     ! Elemental ratios (O/C, H/C, N/C) of the organic species (1 species) in the binary mixture (org+water)
                            & M_org, x_org                  ! Molar mass (g/mol) and mole fraction of the organic species (1 species) in the binary mixture (org+water)
    REAL(rp), PARAMETER :: M_w = 18.01528D0                 ! Molar mass of water (g/mol)
    
    ! Outputs
    REAL(rp)            :: activity_coefficient_water, activity_coefficient_org, &      ! Activity cofficient of water and organic species (1 species) in the binary mixture (org+water)
                            & activity_water, activity_org, &                           ! Activity of water and organic species (1 species) in the binary mixture (org+water)
                            & mass_fraction_water, mass_fraction_org                    ! Mass fraction of water and organic species (1 species) in the binary mixture (org+water)

    
    ! STEP 1: define the organic species (1 species) that is in the binary system (org + water)
    M_org = 200.05D+00
    O2C_org = 0.44D+00
    H2C_org = 1.61D+00  ! set H/C ratio to a negative real value when it is not known
    N2C_org = -1.0D+00  ! set N/C ratio to a negative real value when it is not known
    x_org = 0.75D+00    ! set mole fraction of organic in this example
    
    Group_org = 3
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
      
    ! calculate the activity cofficient of water and organic species:
    BLOCK
        USE BATVBS_MOD  ! Module that contains the main BAT and VBS subroutines
     
        REAL(wrp) :: O2C_org_wrp, H2C_org_wrp, N2C_org_wrp, M_org_wrp, x_org_wrp, &
                     & activity_coefficient_water_wrp, activity_coefficient_org_wrp
        
        ! map precision of REALS: rp -> wrp (defined in PRECISION_MOD)
        M_org_wrp = REAL(M_org, kind = wrp)
        O2C_org_wrp = REAL(O2C_org, kind = wrp)
        H2C_org_wrp = REAL(H2C_org, kind = wrp)
        N2C_org_wrp = REAL(N2C_org, kind = wrp)
        x_org_wrp = REAL(x_org, kind = wrp)        
		
        CALL BAT_SIMULATION(M_org_wrp, O2C_org_wrp, H2C_org_wrp, N2C_org_wrp, Group_org, x_org_wrp, activity_coefficient_water_wrp, &
            & activity_coefficient_org_wrp)
    
        ! map precision of REALS: wrp (defined in PRECISION_MOD) -> rp
        activity_coefficient_water = REAL(activity_coefficient_water_wrp, kind = rp)
        activity_coefficient_org = REAL(activity_coefficient_org_wrp, kind = rp)
    END BLOCK
    
    ! examples of post-processing of BAT calculation output:
    ! calculate the activity of water and organic species (1 species) in the binary mixture (org+water):
    activity_water = activity_coefficient_water * (1.0D+00 - x_org)   ! activity_i = activity_coefficient_i * x_i
    activity_org   = activity_coefficient_org * (x_org)
    
    ! calculate the mass fraction of water and organic species (1 species) in the binary mixture (org+water):
    mass_fraction_water = ((1.0D+00 - x_org) * M_w)/((1.0D+00 - x_org) * M_w + x_org * M_org)     ! mass_fraction_i = (mole_fraction_i x molecular_weight_i)/(mean_molecular_weight)
    mass_fraction_org   =  (x_org * M_org)/((1.0D+00 - x_org) * M_w + x_org * M_org) 
    
    ! example output to screen:
    WRITE(*,'(A)') "Selected predicted properties for given input"
    WRITE(*,'(A,2(ES13.6,1X))') "Activity of water and of organic compound   : ", activity_water, activity_org
    WRITE(*,'(A,2(ES13.6,1X))') "Mass fraction of water and of organic compound: ", mass_fraction_water, mass_fraction_org 
    WRITE(*,'(A)') "Example BAT program completed"
    READ(*,*) ! wait for user action
    
END PROGRAM EXAMPLE_BAT
