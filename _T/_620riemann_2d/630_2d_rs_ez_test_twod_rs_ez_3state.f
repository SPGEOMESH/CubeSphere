














      PROGRAM TEST_TwoD_RS_EZ_3State

      IMPLICIT NONE

      REAL S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1     Ez_1, Ez_2, Ez_3,

     1     Bneta_1, Bneta_2, Bneta_3,

     1     BSSx, BSSy, ESSz,

     1     Bx1, Bx2, Bx3, By1, By2, By3,

     1     Vx_1, Vy_1, Vx_2, Vy_2, Vx_3, Vy_3

      INTEGER i, j, k, ii

      REAL, PARAMETER :: sq3_2 = SQRT ( 3.0) / 2.0

      REAL, DIMENSION ( 3) :: SR, SL, Bx, By, Bneta

      REAL, DIMENSION ( 8) :: V_1, V_2, V_3



      OPEN ( UNIT = 6, FILE = "test_TwoD_RS_EZ_3State.out")



! ----------------------------------------------------------------------
! Initialize the primitive variables "V_1, V_2, V_3" and the

      V_1 ( 1) = 1.0 !Density

      V_1 ( 2) = 0.5 ! Velocity
      V_1 ( 3) = 0.7
      V_1 ( 4) = -0.4

      V_1 ( 5) = 0.8 ! Pressure

      V_1 ( 6) = 1.0 ! Magnetic field
      V_1 ( 7) = 0.3
      V_1 ( 8) = 1.6


! Identical states.
      V_2 = V_1
      V_3 = V_1

! introduce some randomness. All three input states can be different if this is turned on.
#if ( 1 == 1)
      V_2 = V_1 + 0.2
      V_3 = V_1 - 0.4
#endif

      Vx_1 = V_1 ( 2)
      Vy_1 = V_1 ( 3)

      Vx_2 = V_2 ( 2)
      Vy_2 = V_2 ( 3)

      Vx_3 = V_3 ( 2)
      Vy_3 = V_3 ( 3)

! ----------------------------------------------------------------------
! Compute the magnetic field.

      Bx1  = V_1 ( 6) ;     By1 = V_1 ( 7)
      Bx2  = V_2 ( 6) ;     By2 = V_2 ( 7)
      Bx3  = V_3 ( 6) ;     By3 = V_3 ( 7)


      Bneta_1 = - ( By1 + By3) * 0.5                      ! U1, U3

      Bneta_2 = ( sq3_2 * Bx1 + 0.5 * By1
     1          + sq3_2 * Bx2 + 0.5 * By2) * 0.5          ! U1, U2

      Bneta_3 = ( - sq3_2 * Bx2 + 0.5 * By2
     1            - sq3_2 * Bx3 + 0.5 * By3) * 0.5        ! U2, U3

! ----------
! Electric fields.

      Ez_1 = - ( Vx_1 * By1 - Vy_1 * Bx1)
      Ez_2 = - ( Vx_2 * By2 - Vy_2 * Bx2)
      Ez_3 = - ( Vx_3 * By3 - Vy_3 * Bx3)

! ----------
! speeds.

      S_R1 = 1.0 ;       S_L1 = -1.0
      S_R2 = 1.0 ;       S_L2 = -1.0
      S_R3 = 1.0 ;       S_L3 = -1.0

! ------------------------------
! introduce changes to identical states.

! ------------------------------
! test "ST_1, ST_2, ST_3, ST_13, ST_21, ST_32" states.

#if ( 0 == 1)

! #define ST_1

! ----------
#if  ( defined( ST_1 ) )
      S_L1 = 1.0 ;     S_R2 = - 1.0
#elif( defined( ST_2 ) )
      S_L2 = 1.0 ;     S_R3 = - 1.0
#elif( defined( ST_3 ) )
      S_L3 = 1.0 ;     S_R1 = - 1.0
#elif( defined( ST_13 ) )
      S_L3 = 1.0 ;     S_R2 = - 1.0
#elif( defined( ST_21 ) )
      S_L1 = 1.0 ;     S_R3 = - 1.0
#elif( defined( ST_32 ) )
      S_L2 = 1.0 ;     S_R1 = - 1.0
#endif

#endif

! ----------------------------------------------------------------------

      WRITE ( 6, *)"Testing 3-State Two Dimensional Riemann Solver"
      WRITE ( 6, *)" "

! ----------------------------------------

      WRITE ( 6, *)"Input states :"
      WRITE ( 6, *)" "

      DO i = 1, 8

      WRITE ( 6, '(A26, I1, A4)', advance="no")
     1   "V_1, 2, 3, U** ( ", i, ") = "

      WRITE ( 6, '( f8.2, f8.2, f8.2, f10.2 )')
     i   V_1 ( i), V_2 ( i), V_3 ( i)

      END DO

! Magnetic field Components, as seen by 2DRS.
      Bx1 =  ( 2.0 * Bneta_2 + Bneta_1) / SQRT ( 3.0)
      Bx2 =  ( Bneta_2 - Bneta_3) / SQRT ( 3.0)
      Bx3 = -( 2.0 * Bneta_3 + Bneta_1) / SQRT ( 3.0)

      By1 = - Bneta_1
      By2 = Bneta_2 + Bneta_3
      By3 = - Bneta_1


      WRITE ( 6, *)" "
      WRITE ( 6, *)"Magnetic field (as seen by the 2DRS) ",
     1             "and Electric Field:"

      WRITE ( 6, '(A26, I1, A4)', advance="no")" Bx1, By1, Ez_1 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        Bx1, By1, Ez_1

      WRITE ( 6, '(A26, I1, A4)', advance="no")" Bx2, By2, Ez_2 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        Bx2, By2, Ez_2

      WRITE ( 6, '(A26, I1, A4)', advance="no")" Bx3, By3, Ez_3 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        Bx3, By3, Ez_3

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Electrc Fields:"

      WRITE ( 6, '(A26, I1, A4)', advance="no")" Ez_1, Ez_2, Ez_3 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        Ez_1, Ez_2, Ez_3

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Wave speeds:"

      WRITE ( 6, '(A26, I1, A4)', advance="no")" S_R1, S_R2, S_R3 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        S_R1, S_R2, S_R3

      WRITE ( 6, '(A26, I1, A4)', advance="no")" S_L1, S_L2, S_L3 = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')        S_L1, S_L2, S_L3


      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Resolved Fields from Two_D_RS_Ez_3State_HLL:"
      WRITE ( 6, *)" "

      CALL Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"BSSx, BSSy, ESSz = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)')       BSSx, BSSy, ESSz


      END PROGRAM TEST_TwoD_RS_EZ_3State

