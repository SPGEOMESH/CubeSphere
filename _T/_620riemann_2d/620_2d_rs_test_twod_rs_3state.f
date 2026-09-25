




      PROGRAM TEST_TwoD_RS_3State

      IMPLICIT NONE

      INTEGER, PARAMETER :: n_cc_components=8

      REAL n1 ( 3), n2 ( 3), n3 ( 3)

      REAL pi, rhofloor, prsfloor, smallnum, gamma,
     1     isospeed, barotropic_stiffness, gamma_euler,
     1     fourpi, ob_4pi, ob_8pi, E_z, a, c_120, s_120, c_240, s_240,
     1     B0 ( 3)

      REAL, DIMENSION ( n_cc_components) :: U_1, U_2, U_3, F, G,
     1      U_Dbl_Star, F_Dbl_Star, G_Dbl_Star

      REAL, DIMENSION ( n_cc_components) :: V_1, V_2, V_3,
     1                                      W_1, W_2, W_3,
     1                                      F_1, F_2, F_3,
     1                                      G_1, G_2, G_3

      INTEGER i, j, k, n

! ----------------------------------------------------------------------
! Define the constants and initialize values.

      pi = 4.0 * ATAN ( 1.0)

      c_120 = COS ( 120.0 * pi/180.0 )
      s_120 = SIN ( 120.0 * pi/180.0 )

      c_240 = COS ( 240.0 * pi/180.0 )
      s_240 = SIN ( 240.0 * pi/180.0 )

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      gamma = 1.4          ! Polytropic index for the gas.

      isospeed = 1.0              ! Isothermal sound speed.
      barotropic_stiffness = 0.0  ! Stiffens the isothermal EOS if non-zero.
      gamma_euler = 1.1     ! Polytropic index for the gas.

! --------------------------------------------------
! Initialize the conserved variables "U_1, U_2, U_3" and the
! Background magnetic field "B0"

      V_1 ( 1) = 1.0 !Density

      V_1 ( 2) = 0.5 ! Velocity
      V_1 ( 3) = 0.7
      V_1 ( 4) = -0.4

      V_1 ( 5) = 0.8 ! Pressure

      V_1 ( 6) = 1.0 ! Magnetic field
      V_1 ( 7) = 0.3
      V_1 ( 8) = 1.6

! Background magnetic field.
      B0 ( 1)  = 0.0
      B0 ( 2)  = 0.0
      B0 ( 3)  = 0.0

! --------------------------------------------------
! Identical states.

      V_2 = V_1
      V_3 = V_1

! ------------------------------
! introduce changes to identical states.
#if( 1 == 1)

! ------------------------------
! introduce some randomness.
#if 1
      V_2 = V_1 + 0.2
      V_3 = V_1 - 0.4
#endif

! ------------------------------
! test "ST_1, ST_2, ST_3, ST_13, ST_21, ST_32" states.

#define ST_32

! ----------
#if( defined( ST_1 ) )

      V_1 ( 2) = 4.0 * c_240;    V_1 ( 3) = 4.0 * s_240
      V_2 ( 2) = 4.0 * c_240;    V_2 ( 3) = 4.0 * s_240
      V_3 ( 2) = 4.0 * c_240;    V_3 ( 3) = 4.0 * s_240

! ----------
#elif( defined( ST_2 ) )

      V_1 ( 2) = 4.0
      V_2 ( 2) = 4.0
      V_3 ( 2) = 4.0

! ----------
#elif( defined( ST_3 ) )

      V_1 ( 2) = 4.0 * c_120;    V_1 ( 3) = 4.0 * s_120
      V_2 ( 2) = 4.0 * c_120;    V_2 ( 3) = 4.0 * s_120
      V_3 ( 2) = 4.0 * c_120;    V_3 ( 3) = 4.0 * s_120

! ----------
#elif( defined( ST_13 ) )

      V_1 ( 2) = -4.0
      V_2 ( 2) = -4.0
      V_3 ( 2) = -4.0  ;      V_3 = V_1

! ----------
#elif( defined( ST_21 ) )

      V_1 ( 2) = -4.0 * c_120;    V_1 ( 3) = -4.0 * s_120
      V_2 ( 2) = -4.0 * c_120;    V_2 ( 3) = -4.0 * s_120
      V_3 ( 2) = -4.0 * c_120;    V_3 ( 3) = -4.0 * s_120;   V_2 = V_1

! ----------
#elif( defined( ST_32 ) )

      V_1 ( 2) = -4.0 * c_240;    V_1 ( 3) = -4.0 * s_240
      V_2 ( 2) = -4.0 * c_240;    V_2 ( 3) = -4.0 * s_240
      V_3 ( 2) = -4.0 * c_240;    V_3 ( 3) = -4.0 * s_240;   V_3 = V_2

! ----------
#endif

! ------------------------------
! end of introduce changes to identical states.
#endif

! ----------------------------------------------------------------------
! Obtain the conserved vairbles "U_1, U_2, U_3, U_RD" and fluxes,
! "F_1, F_2, F_3, F_RD", "G_1, G_2, G_3, G_RD"
! from the primitive variables "V_1, V_2, V_3, V_RD".

      n = n_cc_components

      CALL MHD_V_TO_UFG ( n, gamma, B0, V_1, U_1, F_1, G_1)
      CALL MHD_V_TO_UFG ( n, gamma, B0, V_2, U_2, F_2, G_2)
      CALL MHD_V_TO_UFG ( n, gamma, B0, V_3, U_3, F_3, G_3)

      W_1 = U_1;  W_2 = U_2;  W_3 = U_3;

! ----------------------------------------------------------------------

!      n1 ( 1) = -0.5;  n1 ( 2) =  SQRT(3.0)/2.0;   n2 ( 3) = 0.0;
!      n2 ( 1) = -0.5;  n2 ( 2) = -SQRT(3.0)/2.0;   n2 ( 3) = 0.0;

      n1 ( 1) =  1.0;  n1 ( 2) =  0.0;             n1 ( 3) = 0.0;
      n2 ( 1) = -0.5;  n2 ( 2) =  SQRT(3.0)/2.0;   n2 ( 3) = 0.0;

      CALL TwoD_RS_3State ( n_cc_components,

     1           pi, rhofloor, prsfloor, smallnum, gamma,
     1           isospeed, barotropic_stiffness, gamma_euler,

     1           n1, n2, n3, B0,

     1           U_1, U_2, U_3,
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

      OPEN ( UNIT = 6, FILE = "test_TwoD_RS_3State.out")

      WRITE ( 6, *)"Testing 3-State Two Dimensional Riemann Solver"
      WRITE ( 6, *)" "

! ----------------------------------------

      WRITE ( 6, *)"Input states and the output U**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "U_1, 2, 3, U** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f10.2 )')
     i   W_1 ( i), W_2 ( i), W_3 ( i), U_Dbl_Star ( i)

      END DO

! ----------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Input states and the output F**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "F_1, 2, 3, F** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f10.2 )')
     i   F_1 ( i), F_2 ( i), F_3 ( i), F_Dbl_Star ( i)

      END DO

! ----------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Input states and the output G**"
      WRITE ( 6, *)" "

      DO i = 1, n_cc_components

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "G_1, 2, 3, G** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f10.2 )')
     i   G_1 ( i), G_2 ( i), G_3 ( i), G_Dbl_Star ( i)

      END DO


      END PROGRAM TEST_TwoD_RS_3State


