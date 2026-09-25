





      SUBROUTINE EVALUATE_TIMESTEP

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           time, dtcur, dtcurp1,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds a new timestep "dtcurp1". We just build the primitive variables
! from the updated conserved variables and use them to build the new
! timestep.

! "U" comes in with just the zone-centered conserved variables and their higher
! order moments. The velocities are defined in a global Cartesian
! space, as are the magnetic fields.

! "V" is the vector of primitive variables that are built from the conserved
! variables.

! "Bxf, Byf, Bzf" are the facial magnetic field components.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

      REAL time, dtcur, dtcurp1

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr

! ----------------------------------------------------------------------

! Variables that are internal to this subroutine.

      INTEGER ixx, iyy, izz, bflag, iprs_flag

      REAL ob_8pi, ob_4pi, a11, a12, a13, a21,
     1     x, y, z, B0_x, B0_y, B0_z, time_use, Btot2,
     1     ent, B0 ( 3), Vpt ( n_cc_components), Upt ( n_cc_components)

! ----------------------------------------------------------------------

!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( CONS_TO_PRIM ) seq

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

#if( CHOOSE_APPL_HELIO == 1)
      dtcurp1 = 1.01 * dtcur
#else
      dtcurp1 = 1.1 * dtcur
#endif

! ------------------------------

! Send the newly initialized variables from CPU to the GPU.

!$acc data copyin ( ob_8pi, ob_4pi)

!$acc update device ( time, dtcur, cfl_coef, dtcurp1)

! ----------------------------------------------------------------------
! We use the primitive variables "V" that are at least second order accurate.
! Also build an estimate of the magnetosonic speed "msonicgr" and use it
! to minimize the next timestep "dtcurp1".


!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, a11, a12, a13, a21,          &
!$OMP&    bflag, time_use, x, y, z, B0_x, B0_y, B0_z, Btot2,            &
!$OMP&    iprs_flag, ent, B0, Upt, Vpt)                                 &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$OMP&    smallnum, ob_4pi, ob_8pi, cfl_coef,                           &
!$OMP&    imask_active_patch, v_scale_length, v_centroid_patch,         &
!$OMP&    U, V, B1_zc, B2_zc, B3_zc,                                    &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, entropygr, msonicgr)       &

!$OMP&            REDUCTION ( MIN: dtcurp1)                             &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$acc&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$acc&    smallnum, ob_4pi, ob_8pi, cfl_coef,                           &
!$acc&   imask_active_patch, v_scale_length, v_centroid_patch,          &
!$acc&   U, V, B1_zc, B2_zc, B3_zc,                                     &
!$acc&    magB_small_cgl, eps_bdry_cgl,                                 &
!$acc&   gamma_eos_gr, gamma_soundspeed_gr, entropygr, msonicgr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, a11, a12, a13, a21,                  &
!$acc&    bflag, time_use, x, y, z, B0_x, B0_y, B0_z, Btot2,            &
!$acc&    iprs_flag, ent, B0, Upt, Vpt)

!$acc&    REDUCTION ( MIN: dtcurp1)

      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) CYCLE

! ------------------------------

          B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)

          x = v_centroid_patch ( ixx, iyy, izz, 1)
          y = v_centroid_patch ( ixx, iyy, izz, 2)
          z = v_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1         ! 1 == magnetic field;  0 == vector potential.
          time_use = 0.0    ! useful for the case of active rotation.

          CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z,
     1                              B0_x, B0_y, B0_z)

          B0 ( 1) = B0_x
          B0 ( 2) = B0_y
          B0 ( 3) = B0_z
#endif

! ------------------------------

          Upt ( :) = U ( ixx, iyy, izz, :, 1)

          ent = entropygr ( ixx, iyy, izz, 1)

          CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)


          U ( ixx, iyy, izz, :, 1) = Upt ( :)
          V ( ixx, iyy, izz, :, 1) = Vpt ( :)

          Btot2 = ( V ( ixx, iyy, izz, 6, 1) + B0 ( 1))**2
     1          + ( V ( ixx, iyy, izz, 7, 1) + B0 ( 2))**2
     1          + ( V ( ixx, iyy, izz, 8, 1) + B0 ( 3))**2

! ------------------------------
#if (IPDE_SYSTEM == 1)
! ------------------------------
#if ( CHOOSE_ADIABATIC == 1)
! Using adiabatic EOS.

          msonicgr ( ixx, iyy, izz) =
     1  gamma_soundspeed_gr ( ixx, iyy, izz) * V ( ixx, iyy, izz, 5, 1)
     1  + ob_4pi * Btot2

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )

#else
! Using isothermal EOS.

          msonicgr ( ixx, iyy, izz) = isospeed**2
     1     + gamma_euler * barotropic_stiffness
     1                   * V ( ixx, iyy, izz, 1, 1)**(gamma_euler - 1.0)
     1     + ob_4pi * Btot2 / V ( ixx, iyy, izz, 1, 1)

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

#endif
! ------------------------------
#elif (IPDE_SYSTEM == 2)
! ------------------------------

          msonicgr ( ixx, iyy, izz) = 1.67 * ( V ( ixx, iyy, izz, 5, 1)
     1                          + 2.0 * V ( ixx, iyy, izz, 9, 1)) / 3.0
     1  + ob_4pi * Btot2

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )

! ------------------------------
#elif (IPDE_SYSTEM == 3)
! ------------------------------

          msonicgr ( ixx, iyy, izz) = 1.67 * ( V ( ixx, iyy, izz, 5, 1)
     1                          + 2.0 * V ( ixx, iyy, izz, 9, 1)) / 3.0
     1                              + 1.67 * V ( ixx, iyy, izz, 10, 1)
     1  + ob_4pi * Btot2

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )

! ------------------------------
#endif
! ------------------------------

          a21 = SQRT ( V ( ixx, iyy, izz, 2, 1)**2
     1               + V ( ixx, iyy, izz, 3, 1)**2
     1               + V ( ixx, iyy, izz, 4, 1)**2)

          a21 = a21 + msonicgr ( ixx, iyy, izz)  ! Maximal signal speed.

          dtcurp1 = AMIN1 ( dtcurp1,
     1    cfl_coef * v_scale_length ( ixx, iyy, izz) / a21 )

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

!$acc update host( dtcurp1 )

! Because "dtcurp1" is a REDUCTION variable in OpenACC, it will be
! automatically copied over to the CPU. No need to do anything further.
! Since, "dtcurp1" is manually moved to the GPU at the beginning program
! there is no harm in updating the CPU value manually here.

! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE EVALUATE_TIMESTEP





