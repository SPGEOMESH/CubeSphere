










      PROGRAM TEST_TwoD_RS_4State

      IMPLICIT NONE

      INTEGER, PARAMETER :: n_cc_components=8

      REAL n1 ( 3), n2 ( 3), n3 ( 3)

      REAL pi, rhofloor, prsfloor, smallnum, gamma,
     1     isospeed, barotropic_stiffness, gamma_euler,
     1     fourpi, ob_4pi, ob_8pi, E_z,
     1     B0 ( 3), a

      REAL, DIMENSION ( n_cc_components) :: U_RU, U_LU, U_LD, U_RD,
     1      U_Dbl_Star, F_Dbl_Star, G_Dbl_Star, F, G

      REAL, DIMENSION ( n_cc_components) :: V_RU, V_LU, V_LD, V_RD,
     1                                      W_RU, W_LU, W_LD, W_RD,
     1                                      F_RU, F_LU, F_LD, F_RD,
     1                                      G_RU, G_LU, G_LD, G_RD

      INTEGER i, j, k, n

! ----------------------------------------------------------------------
! Define the constants and initialize values.

      pi = 4.0 * ATAN ( 1.0)

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      gamma = 1.4          ! Polytropic index for the gas.

      isospeed = 1.0              ! Isothermal sound speed.
      barotropic_stiffness = 0.0  ! Stiffens the isothermal EOS if non-zero.
      gamma_euler = 1.1     ! Polytropic index for the gas.

! --------------------------------------------------
! Initialize the conserved variables "U_RU, U_LU, U_LD, U_RD" and the
! Background magnetic field "B0"

      V_RU ( 1) = 1.0 !Density

      V_RU ( 2) = 0.5 ! Velocity
      V_RU ( 3) = 0.7
      V_RU ( 4) = -0.4

      V_RU ( 5) = 0.8 ! Pressure

      V_RU ( 6) = 1.0 ! Magnetic field
      V_RU ( 7) = 0.3
      V_RU ( 8) = 1.6

! Background magnetic field.
      B0 ( 1)  = 0.0
      B0 ( 2)  = 0.0
      B0 ( 3)  = 0.0

! --------------------------------------------------
! Identical states.

      V_LU = V_RU
      V_LD = V_RU
      V_RD = V_RU

! ------------------------------
! introduce changes to identical states.
#if( 1 == 1)

! ------------------------------
! introduce some randomness.
#if 1
      V_LU = V_RU + 0.2
      V_LD = V_RU - 0.4
      V_RD = V_RU + 1.0
#endif

! ------------------------------
! test DOWN, UP, RIGHT, LEFT, RU, LU, LD, RD states.

#define RD

! ----------
#if( defined( DOWN ) )

      V_LD = V_RD
      V_RU ( 3) = 4.0
      V_LU ( 3) = 4.0
      V_LD ( 3) = 4.0
      V_RD ( 3) = 4.0

! ----------
#elif( defined( UP ) )

      V_LU = V_RU
      V_RU ( 3) = -4.0
      V_LU ( 3) = -4.0
      V_LD ( 3) = -4.0
      V_RD ( 3) = -4.0

! ----------
#elif( defined( LEFT ) )

      V_LU = V_LD
      V_RU ( 2) = 4.0
      V_LU ( 2) = 4.0
      V_LD ( 2) = 4.0
      V_RD ( 2) = 4.0

! ----------
#elif( defined( RIGHT ) )

      V_RU = V_RD
      V_RU ( 2) = -4.0
      V_LU ( 2) = -4.0
      V_LD ( 2) = -4.0
      V_RD ( 2) = -4.0

! ----------
#elif( defined( RU ) )

      V_RU ( 2) = -4.0;   V_RU ( 3) = -4.0
      V_LU ( 2) = -4.0;   V_LU ( 3) = -4.0
      V_LD ( 2) = -4.0;   V_LD ( 3) = -4.0
      V_RD ( 2) = -4.0;   V_RD ( 3) = -4.0

! ----------
#elif( defined( LU ) )

      V_RU ( 2) =  4.0;   V_RU ( 3) = -4.0
      V_LU ( 2) =  4.0;   V_LU ( 3) = -4.0
      V_LD ( 2) =  4.0;   V_LD ( 3) = -4.0
      V_RD ( 2) =  4.0;   V_RD ( 3) = -4.0

! ----------
#elif( defined( LD ) )

      V_RU ( 2) =  4.0;   V_RU ( 3) =  4.0
      V_LU ( 2) =  4.0;   V_LU ( 3) =  4.0
      V_LD ( 2) =  4.0;   V_LD ( 3) =  4.0
      V_RD ( 2) =  4.0;   V_RD ( 3) =  4.0

! ----------
#elif( defined( RD ) )

      V_RU ( 2) = -4.0;   V_RU ( 3) =  4.0
      V_LU ( 2) = -4.0;   V_LU ( 3) =  4.0
      V_LD ( 2) = -4.0;   V_LD ( 3) =  4.0
      V_RD ( 2) = -4.0;   V_RD ( 3) =  4.0

! ----------
#endif

! ------------------------------
! end of introduce changes to identical states.
#endif

! ----------------------------------------------------------------------
! Obtain the conserved vairbles "U_RU, U_LU, U_LD, U_RD" and fluxes,
! "F_RU, F_LU, F_LD, F_RD", "G_RU, G_LU, G_LD, G_RD"
! from the primitive variables "V_RU, V_LU, V_LD, V_RD".

      n = n_cc_components

      CALL MHD_V_TO_UFG ( n, gamma, B0, V_RU, U_RU, F_RU, G_RU)
      CALL MHD_V_TO_UFG ( n, gamma, B0, V_LU, U_LU, F_LU, G_LU)
      CALL MHD_V_TO_UFG ( n, gamma, B0, V_LD, U_LD, F_LD, G_LD)
      CALL MHD_V_TO_UFG ( n, gamma, B0, V_RD, U_RD, F_RD, G_RD)

      W_RU = U_RU;  W_LU = U_LU;  W_LD = U_LD;  W_RD = U_RD;

! ----------------------------------------------------------------------

!     n1 ( 1) = 1.0;  n1 ( 2) =  0.0;   n1 ( 3) = 0.0;
!     n2 ( 1) = -0.6;  n2 ( 2) =  1.0;   n2 ( 3) = 0.0;

      n1 ( 1) =  1.0;  n1 ( 2) =   0.3;    n1 ( 3) = 0.0;
      n2 ( 1) = -0.5;  n2 ( 2) =   1.0;    n2 ( 3) = 0.0;

      CALL TwoD_RS_4State ( n_cc_components,

     1           pi, rhofloor, prsfloor, smallnum, gamma,
     1           isospeed, barotropic_stiffness, gamma_euler,

     1           n1, n2, n3, B0,

     1           U_RU, U_LU, U_LD, U_RD,
     1           U_Dbl_Star, F_Dbl_Star, G_Dbl_Star, E_z )

! Rotate from global x-coordinate to the "n1" axis. This is only needed
! to check the correctness of the 2D RS.
      a = ATAN2 ( n1 ( 2), n1 ( 1))

      F = F_Dbl_Star
      G = G_Dbl_Star

      F_Dbl_Star = F * COS ( a) - G * SIN ( a)
      G_Dbl_Star = F * SIN ( a) + G * COS ( a)

! ----------------------------------------------------------------------
! Testing for orthonormality.
! Declare the 3x3 arrays "t_n, t_np, t_nxnp" in the subrotine
! and paste the following lines inside the subroutine, for testing.
! It is expected to write out an identity matrix.

#if( 1 == 0)

      t_n ( 1, 1) = n1( 1); t_n ( 1, 2) = n1( 2); t_n ( 1, 3) = n1( 3)
      t_n ( 2, 1) = n12( 1);t_n ( 2, 2) = n12( 2);t_n ( 2, 3) = n12( 3)
      t_n ( 3, 1) = n3( 1); t_n ( 3, 2) = n3( 2); t_n ( 3, 3) = n3( 3)

      t_np( 1, 1) = np1( 1); t_np( 1, 2) = np1( 2); t_np( 1, 3) = np1(3)
      t_np( 2, 1) = np2( 1); t_np( 2, 2) = np2( 2); t_np( 2, 3) = np2(3)
      t_np( 3, 1) = np3( 1); t_np( 3, 2) = np3( 2); t_np( 3, 3) = np3(3)

      t_nxnp = MATMUL ( t_n, t_np)

      WRITE( *,*)"t_nxnp",t_nxnp

#endif

! ----------------------------------------------------------------------

      OPEN ( UNIT = 6, FILE = "test_TwoD_RS_4State.out")

      WRITE ( 6, *)"Testing 4-State Two Dimensional Riemann Solver"
      WRITE ( 6, *)" "

! ----------------------------------------

      WRITE ( 6, *)"Input states and the output U**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "U_RU, LU, LD, RD, U** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f8.2, f10.2 )')
     i   W_RU ( i), W_LU ( i), W_LD ( i), W_RD ( i), U_Dbl_Star ( i)

      END DO

! ----------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Input states and the output F**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "F_RU, LU, LD, RD, F** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f8.2, f10.2 )')
     i   F_RU ( i), F_LU ( i), F_LD ( i), F_RD ( i), F_Dbl_Star ( i)

      END DO

! ----------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Input states and the output G**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "G_RU, LU, LD, RD, G** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f8.2, f10.2 )')
     i   G_RU ( i), G_LU ( i), G_LD ( i), G_RD ( i), G_Dbl_Star ( i)

      END DO

      END PROGRAM TEST_TwoD_RS_4State


