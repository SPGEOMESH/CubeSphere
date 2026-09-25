














      PROGRAM TEST_TwoD_RS_EZ_4State

      IMPLICIT NONE

      REAL sin_t, cos_t, S_R, S_L, S_U, S_D,
     1     Ez_RU, Ez_LU, Ez_LD, Ez_RD,
     1     Bneta_U, Bneta_D, By_R, By_L,
     1     BSSx, BSSy, ESSz,

     1     Bx_RU, By_RU, Bx_LU, By_LU, Bx_LD, By_LD, Bx_RD, By_RD,
     1     Vx_RU, Vy_RU, Vx_LU, Vy_LU, Vx_LD, Vy_LD, Vx_RD, Vy_RD

      INTEGER i, j, k, ii

      OPEN ( unit = 6, file = "test_TwoD_RS_EZ_4State.out")



      S_R = 1.0 ; S_L = - 1.0 ; S_U = 1.0 ; S_D = -1.0

#if ( 1 == 1)

      sin_t = 1.0 ; cos_t = SQRT ( 1.0 - sin_t**2)

      Bneta_U = 1.0 ; Bneta_D = 1.0 ; By_R = 0.27 ; By_L = 0.27

      Vx_RU = 0.5 ; Vy_RU = - 0.3 ; Vx_LU = 0.5 ; Vy_LU = - 0.3
      Vx_LD = 0.5 ; Vy_LD = - 0.3 ; Vx_RD = 0.5 ; Vy_RD = - 0.3

#else

      sin_t = 0.7 ; cos_t = SQRT ( 1.0 - sin_t**2)

      Bneta_U = 1.2 ; Bneta_D = 0.6 ; By_R = 0.2 ; By_L = 0.7

      Vx_RU = 0.2 ; Vy_RU = - 0.1 ; Vx_LU = 0.4 ; Vy_LU = - 0.3
      Vx_LD = 0.5 ; Vy_LD = - 0.4 ; Vx_RD = 0.6 ; Vy_RD = - 0.2

#endif

      By_RU = By_R  ;  Bx_RU = ( Bneta_U + cos_t * By_RU) / sin_t
      By_LU = By_L  ;  Bx_LU = ( Bneta_U + cos_t * By_LU) / sin_t
      By_RD = By_R  ;  Bx_RD = ( Bneta_D + cos_t * By_RD) / sin_t
      By_LD = By_L  ;  Bx_LD = ( Bneta_D + cos_t * By_LD) / sin_t

      Ez_RU = - ( Vx_RU * By_RU - Vy_RU * Bx_RU )
      Ez_LU = - ( Vx_LU * By_LU - Vy_LU * Bx_LU )
      Ez_LD = - ( Vx_LD * By_LD - Vy_LD * Bx_LD )
      Ez_RD = - ( Vx_RD * By_RD - Vy_RD * Bx_RD )

999   FORMAT ( 1pe12.4, 2x, 1pe12.4, 2x, 1pe12.4, 2x, 1pe12.4 )
998   FORMAT ( 1pe12.4, 2x, 1pe12.4, 2x, 1pe12.4 )

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sin_t, cos_t = ", sin_t, cos_t

      WRITE ( 6, *)" "
      WRITE ( 6, *)"S_R, S_L, S_U, S_D = "
      WRITE ( 6, 999) S_R, S_L, S_U, S_D

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Bneta_U, Bneta_D, By_R, By_L = "
      WRITE ( 6, 999) Bneta_U, Bneta_D, By_R, By_L

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Bx_RU, By_RU, Bx_LU, By_LU = "
      WRITE ( 6, 999) Bx_RU, By_RU, Bx_LU, By_LU

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Bx_RD, By_RD, Bx_LD, By_LD = "
      WRITE ( 6, 999) Bx_RD, By_RD, Bx_LD, By_LD

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Ez_RU, Ez_LU, Ez_LD, Ez_RD = "
      WRITE ( 6, 999) Ez_RU, Ez_LU, Ez_LD, Ez_RD

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Calling Two_D_RS_Ez_4State_HLL"

      CALL Two_D_RS_Ez_4State_HLL

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"BSSx, BSSy, ESSz = "
      WRITE ( 6, 998) BSSx, BSSy, ESSz


      END PROGRAM TEST_TwoD_RS_EZ_4State



