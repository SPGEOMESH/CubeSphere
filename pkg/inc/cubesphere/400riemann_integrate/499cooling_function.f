





! All stuff related to cooling function
! that is wrapped in USE_COOLING_FUNCTION.
! That is all I know.


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( USE_COOLING_FUNCTION == 1)

      SUBROUTINE OPERATOR_SPLIT_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1           U, entropygr, Bxf, Byf, Bzf, flag_Tfloor)

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine does an operator splitting on the cooling source term.
! The internal energy and hence the pressure is evolved, based on the
! amount of cooling provided by the cooling function.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

      REAL time, dtcur

! ------------------------------
! VARIABLES THAT GENERATE THE MESH (set in EXTRACT_VERT_COORDS_PATCH):-

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------
! VOLUMETRIC VARIABLES (set in COORDS_TO_WT_QUADPTS_VOLUME_ISO2):-

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------
! FLUID VARIABLES:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

! ------------------------------
! Extra variable for silo output; tempfloor flag.

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost) ::
     1   flag_Tfloor

! ----------------------------------------------------------------------
! Variables past this point are local to this subroutine.

      REAL rho, prs0, prs1, dt_cool,
     1     vx, vy, vz, Bx, By, Bz, a1, a2, ent, ob_8pi,
     1     is_tfloor, rad

      INTEGER ixx, iyy, izz, j, iexpand

! ----------------------------------------------------------------------

!$acc routine( COOLING_PTWISE ) seq

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)

!$acc data copyin ( ob_8pi)

!$acc kernels
      flag_Tfloor ( :, :, :) = 0.0
!$acc end kernels

! ----------------------------------------------------------------------
! Operator splitting source term: cooling function.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, rho, prs0, prs1, dt_cool,    &
!$OMP&      vx, vy, vz, Bx, By, Bz, a1, a2, ent, is_tfloor, rad)        &
!$OMP&          SHARED ( nr_zones, ntp_zones, NGhost, j, dtcur,         &
!$OMP&      gamma, rhofloor, prsfloor, gasconst, molewt, smallnum, pi,  &
!$OMP&      isospeed, gamma_euler, barotropic_stiffness, tempfloor,     &
!$OMP&      v_centroid_patch,                                           &
!$OMP&      ob_8pi, U, entropygr, flag_Tfloor)                          &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present ( dtcur, ob_8pi,               &
!$acc&      gamma, rhofloor, prsfloor, gasconst, molewt, smallnum, pi,  &
!$acc&      isospeed, gamma_euler, barotropic_stiffness, tempfloor,     &
!$acc&      v_centroid_patch,                                           &
!$acc&      U, entropygr, flag_Tfloor)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, rho, prs0, prs1, dt_cool,            &
!$acc&      vx, vy, vz, Bx, By, Bz, a1, a2, ent, is_tfloor, rad)

      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

! --------------------------------------------------

          rad = SQRT ( v_centroid_patch ( ixx, iyy, izz, 1)**2
     1               + v_centroid_patch ( ixx, iyy, izz, 2)**2
     1               + v_centroid_patch ( ixx, iyy, izz, 3)**2)

          rho = MAX ( U ( ixx, iyy, izz, 1, 1), rhofloor)

          a1 = rho / ( rho**2 + rhofloor**2)

          vx = U ( ixx, iyy, izz, 2, 1) * a1
          vy = U ( ixx, iyy, izz, 3, 1) * a1
          vz = U ( ixx, iyy, izz, 4, 1) * a1

          Bx = U ( ixx, iyy, izz, 6, 1)
          By = U ( ixx, iyy, izz, 7, 1)
          Bz = U ( ixx, iyy, izz, 8, 1)

          a1 = ( vx**2 + vy**2 + vz**2) * 0.5 * rho
          a2 = ( Bx**2 + By**2 + Bz**2) * ob_8pi

          ent = entropygr ( ixx, iyy, izz, 1)

! ----------
! check if the floor temperature is used to reset.

          is_tfloor = 0.0

! ----------
! Cooling is done only for the adiabatic case.

#if ( CHOOSE_ADIABATIC == 1)

          prs0 = ( gamma - 1.0 )
     1         * ( U ( ixx, iyy, izz, 5, 1) - a1 - a2)

          IF ( prs0 <= prsfloor) THEN
#if ( PRESSURE_POSITIVITY >= 1)
            prs0 = MAX ( prsfloor, ent * rho**( gamma - 1.0))
#else
            prs0 = prsfloor
#endif
          END IF

! ----------

          CALL COOLING_PTWISE

     1   ( pi, rhofloor, prsfloor, tempfloor, smallnum,
     1     gamma, molewt, gasconst,
     1     isospeed, barotropic_stiffness, gamma_euler, rad,

     1     dtcur, rho, prs0, prs1, dt_cool, is_tfloor)

! ----------
! Update the energy with the evolved "prs1".

          U ( ixx, iyy, izz, 5, 1) = a1 + a2 + ( prs1 / ( gamma - 1.0))

          flag_Tfloor ( ixx, iyy, izz) = is_tfloor

#endif
! #if ( CHOOSE_ADIABATIC == 1)
! --------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE OPERATOR_SPLIT_SOURCE





      SUBROUTINE COOLING_PTWISE

     1   ( pi, rhofloor, prsfloor, tempfloor, smallnum,
     1     gamma, molewt, gasconst,
     1     isospeed, barotropic_stiffness, gamma_euler, rad,

     1     dtcur, rho, prs0, prs1, dt_cool, is_tfloor)

! ----------------------------------------------------------------------
! Written By: Sethupathy S.
!
! The cooling is applied to the energy equation using the operator-splitting.
!
! In this subroutine, the ODE
!               "dy / dt = -Lm * rho**2" equation is solved.
!
! "y" is the internal energy; "Lm" is the scaled cooling function;
! "Lm * rho**2" is the cooling power.
!
! "prs0" is the initial pressure.
! "prs1" is the updated pressure.
!
! If the temperature is less than the star temperature; then there is no cooling.
! This is enforced in the subroutine GET_COOLING_FUNCTION_VALUE ( ).
!
! ----------------------------------------------------------------------
!$acc routine seq

      IMPLICIT NONE

      REAL pi, rhofloor, prsfloor, tempfloor, smallnum,
     1     gamma, molewt, gasconst,
     1     isospeed, barotropic_stiffness, gamma_euler, rad,

     1     dtcur, rho, prs0, prs1, dt_cool, is_tfloor

! ----------------------------------------------------------------------
! Variables past this point are local to this subroutine.

      INTEGER i, j, ii

      REAL Tfloor, T, Lm, dLmdT, prs0_fl


! Cash-Karp parameters for the adaptive stepsize, taken from Chapter 16.2
! of Numerical recipes book.

      REAL, PARAMETER ::
     1      c1 =  37.0/378.0,  c2 = 0.0,  c3 = 250.0/ 621.0,
     1      c4 = 125.0/594.0,  c5 = 0.0,  c6 = 512.0/1771.0,

     1      cs1 =  2825.0/27648.0,  cs2 = 0.0,  cs3 = 18575.0/48384.0,
     1      cs4 = 13525.0/55296.0,              cs5 =   277.0/14336.0,
     1      cs6 =     1.0/4.0,

     1      b21 = 1.0/ 5.0,
     1      b31 = 3.0/40.0,  b32 =  9.0/40.0,
     1      b41 = 3.0/10.0,  b42 = -9.0/10.0,  b43 = 6.0/5.0,

     1      b51 = -11.0/54.0,  b52 = 5.0/2.0,  b53 = -70.0/27.0,
     1      b54 =  35.0/27.0,

     1      b61 = 1631.0/55296.0,  b62 = 175.0/512.0,
     1      b63 =  575.0/13824.0,  b64 = 44275.0/110592.0,
     1      b65 =  253.0/ 4096.0

! --------------------

      REAL k1, k2, k3, k4, k5, k6, f, y0, y1, ys1, prs, Uph,

     1     ta, tf, dt, eng_to_T, T0, T1, scl_err

! ----------------------------------------------------------------------

!$acc routine( GET_COOLING_FUNCTION_VALUE ) seq

! ----------------------------------------------------------------------

      prs0_fl = 0.0        ! rho * gasconst * tempfloor / molewt

! start and final time for the adaptive stepping.
      ta = 0.0
      tf = dtcur
      dt = dtcur

      Tfloor = tempfloor

      eng_to_T = ( gamma - 1.0) * molewt / ( gasconst * rho)

      y0 = prs0 / ( gamma - 1.0)         ! Internal energy.

      T0 = y0 * eng_to_T                 ! Temperature.

      prs1 = prs0

! --------------------------------------------------
! reach "tf" in 1000 or less iterations.

      DO i = 1, 1000

! ----------------------------------------
! count the iterations.

       ii = i

! ------------------------------

        y1 = y0

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k1 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

        y1 = y0 + b21*k1

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k2 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

        y1 = y0 + b31*k1 + b32*k2

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k3 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

        y1 = y0 + b41*k1 + b42*k2 + b43*k3

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k4 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

        y1 = y0 + b51*k1 + b52*k2 + b53*k3 + b54*k4

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k5 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

        y1 = y0 + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5

! ----------

        T   = y1 * eng_to_T        ! Temperature.
        prs = y1 * ( gamma - 1.0)  ! Pressure.
        prs = MAX ( prs, prs0_fl)

        CALL GET_COOLING_FUNCTION_VALUE( Tfloor, T, Lm, dLmdT, rad, Uph)

        k6 = dt * ( - Lm * rho**2 - Uph * ( prs - prs0_fl))

! --------------------

        y1 = y0 + c1*k1 + c2*k2 + c3*k3
     1          + c4*k4 + c5*k5 + c6*k6       ! 5th order RK update.

        ys1 = y0 + cs1*k1 + cs2*k2 + cs3*k3
     1           + cs4*k4 + cs5*k5 + cs6*k6   ! 4th order update.

! --------------------
! scaled error.

        scl_err = ABS ( y1 - ys1) / ( ABS ( y0) + ABS ( k1))

        scl_err = MAX ( scl_err / 1.0e-5, 0.0)

! ----------------------------------------

        IF ( scl_err < 1.0) THEN

! --------------------
! if the scaled error is less than 1.0 then the "y1" is considered accurate.

          y0   = y1                            ! update internal energy.
          prs1 = y1 * ( gamma - 1.0)           ! update pressure.
          ta   = ta + dt                       ! update time.

          IF ( ABS ( ta / tf - 1.0 ) < 1.0e-9) EXIT   ! exit at "ta == tf".

          dt = MIN ( dt * 2.0,                 ! increase "dt".
     1               dt * scl_err**( - 0.2),
     1               tf - ta)

! --------------------

        ELSE

! --------------------
! calculated "y1" is not accurate; reduce "dt" and calculate "y1" again.

          dt = MAX ( dt * 0.05,
     1               dt * 0.8 * scl_err**( - 0.25))

! --------------------

        END IF

! ----------------------------------------

      END DO  ! End of "DO i = 1, 1000" loop.

! ----------------------------------------------------------------------
! Enforce temperature floor.

      T1   = y1 * eng_to_T                 ! Temperature
      prs1 = y1 * ( gamma - 1.0)           ! Pressure

!     IF ( ( T1 .LE. tempfloor) .AND. ( T0 .GE. tempfloor)) THEN
      IF ( ( T1 .LE. tempfloor)) THEN

        is_tfloor = 1.0

        y1   = tempfloor / eng_to_T           ! internal energy

        prs1 = y1 * ( gamma - 1.0)

      END IF

! ----------------------------------------------------------------------

#if ( 0 == 1)

      WRITE ( *, *)"No. of steps = ", ii

#else

      IF ( ii .GE. 1000) THEN
        WRITE ( *, *)"No. of steps in COOLING_PTWISE ( )= ", ii
        STOP
      END IF

#endif

! ----------------------------------------------------------------------

      END SUBROUTINE COOLING_PTWISE

#endif

