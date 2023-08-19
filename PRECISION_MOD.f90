!****************************************************************************************
!*   MODULE PRECISION_MOD                                                               *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Module that defines the precision of REAL/INTEGER data type, the length of         *
!*   CHAR type and a small number used in conditional statements.                       *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes: 06-05-2021                                                      *
!****************************************************************************************

MODULE PRECISION_MOD

    USE, INTRINSIC :: IEEE_ARITHMETIC ! Intrinsic module used to specify the minimum number of decimal digits for reals

    IMPLICIT NONE

    INTEGER(4), PARAMETER, PRIVATE :: number_of_digits = 15                           ! Desired minimum number of precision for reals
    INTEGER(4), PARAMETER, PUBLIC  :: wrp = IEEE_SELECTED_REAL_KIND(number_of_digits) ! Working precision level of reals
    INTEGER(4), PARAMETER, PUBLIC  :: wip = 4                                         ! Working precision level of integers
    INTEGER(4), PARAMETER, PUBLIC  :: wlp = 4                                         ! Working precision level of logical values
    INTEGER(4), PARAMETER, PUBLIC  :: clen = 45                                       ! Length of character type

    REAL(wrp), PARAMETER, PUBLIC :: tinynumber = EPSILON(1.0_wrp)                     ! Used to test == 0.0_wrp,
                                                                                      ! < 0.0_wrp, > 0.0_wrp, etc.

END MODULE PRECISION_MOD
