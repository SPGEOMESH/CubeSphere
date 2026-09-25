




      PROGRAM TEST_COOLING_PTWISE

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

#include <cubesphere/000com/cool.com>

      INTEGER, PARAMETER :: N = 100

      REAL G_iprob, M_iprob, V_iprob, B_iprob, cs_iprob,
     1     om_iprob, angle_iprob, alph_iprob, Tfloor_iprob,

     1     rho, prs0, prs1, dtcur, dt_cool, is_tfloor,

     1     pow, x, y, z, V_CAK_radial ( 3), rhoc_o, T0, Lm, dLmdT,
     1     rad, Uph,

     1     rmin, rmax

      REAL radius ( 0: N)

      INTEGER i, j, ii

! ----------------------------------------------------------------------

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

! Initialize "rmin, rmax" to reasonable values.
      rmin = 1.0    ! 10.0
      rmax = 15.0   ! 15.23

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 0.62   ! Set up mole wt. of the gas.

      gamma_euler = 5.0 / 3.0    ! Polytropic index for the gas.
      gasconst_euler = 1.0       ! Set up scaled gas constant.
      molewt_euler = 0.62        ! Set up mole wt. of the gas.

      isospeed = cs_iprob         ! Isothermal sound speed; only for isomhd problems.
      barotropic_stiffness = 0.0  ! Stiffens the isothermal EOS if non-zero.

      tempfloor = 1.0e-5          ! used only when the cooling fn is called.

! ----------------------------------------------------------------------

! For CAK winds problem:-

      CALL GET_CAK_PARAMETERS( omega_passive, isospeed,
     1                         dipole_B, alph_cak,
     1                         Mstar, newtonsc, tempfloor,
     1                         iflag_npoles, multipole_Bmag)

      WRITE( 6, *)" "
      WRITE( 6, *)"omega_passive = ", omega_passive
      WRITE( 6, *)"dipole_B = ", dipole_B
      WRITE( 6, *)"isospeed, Mstar, newtonsc, tempfloor = ",
     1             isospeed, Mstar, newtonsc, tempfloor
      WRITE( 6, *)"iflag_npoles   = ", iflag_npoles
      WRITE( 6, *)"multipole_Bmag = ", multipole_Bmag
      WRITE( 6, *)" "


! Fill in the cooling table. different cooling functions can be read in by making
! small modifications to this subroutine.
      CALL FILL_CAK_COOLTABLE ( )

! ----------------------------------------------------------------------

      dtcur = 1.0e-6

      DO i = 0, N

        pow = REAL ( i) / REAL ( N)
        radius ( i) = rmin * ( rmax / rmin)**pow

      END DO

! ----------------------------------------------------------------------

      DO i = 0, N

! --------------------------------------------------

      x = radius ( i)
      y = 0.0
      z = 0.0

      rad = radius ( i)

      CALL CAK_INIT_1D ( x, y, z, V_CAK_radial, rhoc_o)

      rho  = V_CAK_radial ( 1)
      prs0 = V_CAK_radial ( 3)

! ----------

!     T0 = prs0 * molewt / ( gasconst * rho)     ! Temperature

      IF ( MOD ( i, 5) == 0) T0 = 2.0
      IF ( MOD ( i, 5) == 1) T0 = 10.0
      IF ( MOD ( i, 5) == 2) T0 = 20.0
      IF ( MOD ( i, 5) == 3) T0 = 50.0
      IF ( MOD ( i, 5) == 4) T0 = 20000.0

      prs0 = T0 * ( gasconst * rho) / molewt

      CALL GET_COOLING_FUNCTION_VALUE ( tempfloor, T0, Lm, dLmdT,
     1                                  rad, Uph)

! ----------
! check the usage of temperature floor.

      is_tfloor = 0.0

! ----------

      CALL COOLING_PTWISE

     1   ( pi, rhofloor, prsfloor, tempfloor, smallnum,
     1     gamma, molewt, gasconst,
     1     isospeed, barotropic_stiffness, gamma_euler, rad,

     1     dtcur, rho, prs0, prs1, dt_cool, is_tfloor)

! ----------

!     WRITE( 6, *)"i, R, T0, Lm, prs0, prs1 = ",
!    1   i, radius ( i), T0, Lm, prs0, prs1

      WRITE( 6, *)"i, R, T0, Lm, T1 = ",
     1   i, radius ( i), T0, Lm, prs1 * molewt / ( gasconst * rho)

! --------------------------------------------------

      END DO   ! End of "DO i = 1, N" loop.

! ----------------------------------------------------------------------

      END PROGRAM TEST_COOLING_PTWISE

