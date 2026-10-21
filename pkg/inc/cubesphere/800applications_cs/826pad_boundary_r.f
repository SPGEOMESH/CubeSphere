









      SUBROUTINE PAD_BOUNDARY_R

     1    ( nr_zones, ntp_zones, NGhost,
     1      n_cc_components, n_cc_modes,
     1      n_fc_modes,

     1      vert_coords_patch, imask_active_patch,

     1      xf_area_patch, yf_area_patch, zf_area_patch,
     1      v_volume_patch,

     1      v_centroid_patch,

     1      U, entropygr, Bxf, Byf, Bzf,

     1      bcarr_r)

! ----------------------------------------------------------------------
! Written by : Sethupathy S.

! Pads the physical inner or outer boundaries in the radial direction.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

! Arrays describing the geometry of the mesh.

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      REAL xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! Arrays representing the MHD variables in volume and on faces.

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

      INTEGER, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost, 2) :: bcarr_r

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, i, j, ilhs, irhs, oned_begin, oned_end

      REAL ob_8pi

! ----------------------------------------------------------------------
! Application specific variables.

      REAL xu, yu, zu, rlen, rho_vr           ! used in CAK simulations.

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * 4.0 * ATAN ( 1.0))

!$acc data copyin( ob_8pi )

! ----------------------------------------------------------------------
#if ( CHOOSE_APPL_HELIO == 1)

! This boundary condition has only been necessary for helosphere problem.

! For the case with the inflow/outflow velocity at the outer boudnary, assign the
! correct "bcarr_r" for each zone at the outer boundary.

! In other words, for "bcarr_r ( ixx, iyy, 2) == 3" check each zone and
! assign,
!      bcarr_r ( ixx, iyy, 2) =  2 ;  If ( u.r >= 0)
!  or, bcarr_r ( ixx, iyy, 2) = -1 ;  if ( u.r <  0)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, xu, yu, zu, rlen, irhs)           &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    imask_active_patch, inflow_uvec, bcarr_r, v_centroid_patch)
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    imask_active_patch, inflow_uvec, bcarr_r, v_centroid_patch)
!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, xu, yu, zu, rlen, irhs)

      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      IF ( imask_active_patch ( ixx, iyy, 2) .EQ. 1) THEN

      IF ( bcarr_r ( ixx, iyy, 2) == 3 ) THEN

! --------------------

      irhs = nr_zones

      rlen = SQRT ( v_centroid_patch ( ixx, iyy, irhs, 1)**2
     1            + v_centroid_patch ( ixx, iyy, irhs, 2)**2
     1            + v_centroid_patch ( ixx, iyy, irhs, 3)**2)

      xu = v_centroid_patch ( ixx, iyy, irhs, 1) / rlen
      yu = v_centroid_patch ( ixx, iyy, irhs, 2) / rlen
      zu = v_centroid_patch ( ixx, iyy, irhs, 3) / rlen

      IF ( ( xu * inflow_uvec ( 1)
     1     + yu * inflow_uvec ( 2)
     1     + zu * inflow_uvec ( 3)) .GE. 0.0) THEN

        bcarr_r ( ixx, iyy, 2) = 2

      ELSE

        bcarr_r ( ixx, iyy, 2) = -1   ! Outflow NOT assigned.

      END IF  ! End of "IF ( ( xu * inflow_uvec ( 1) ...".

! --------------------

      END IF   ! End of "IF ( bcarr_r ( ixx, iyy, 2) == 3 ) THEN".
      END IF   ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO
      END DO

!$acc end parallel
!$OMP END PARALLEL DO

#endif
! End of "( CHOOSE_APPL_HELIO == 1)"

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ilhs, irhs)                  &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$OMP&    v_volume_patch, v_centroid_patch,                             &
!$OMP&    imask_active_patch, bcarr_r, U, Bxf, Byf, Bzf)                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$acc&    v_volume_patch, v_centroid_patch,                             &
!$acc&   imask_active_patch, bcarr_r, U, Bxf, Byf, Bzf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ilhs, irhs)

      DO izz = 1, NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

! ------------------------------

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ------------------------------
! Outflow/Continuitive at lower r-boundary.

      IF ( bcarr_r ( ixx, iyy, 1) == 2 ) THEN

! --------------------

      ilhs = 1 - izz
      irhs = 1

! ----------
! Inner boundary is not pad-boundary for practical problems.
! So, this is practically never used.

! ----------
#if( 0 == 1)
! ----------

      U ( ixx, iyy, ilhs, :, 1) = U ( ixx, iyy, irhs, :, 1)
     1             * v_volume_patch ( ixx, iyy, irhs)
     1             / v_volume_patch ( ixx, iyy, ilhs)

      Bzf ( ixx, iyy, ilhs-1, 1) = Bzf ( ixx, iyy, irhs-1, 1)
     1                 * zf_area_patch ( ixx, iyy, irhs-1)
     1                 / zf_area_patch ( ixx, iyy, ilhs-1)

      Bxf ( ixx, iyy, ilhs, 1) = Bxf ( ixx, iyy, irhs, 1)
     1               * xf_area_patch ( ixx, iyy, irhs)
     1               / xf_area_patch ( ixx, iyy, ilhs)

      Byf ( ixx, iyy, ilhs, 1) = Byf ( ixx, iyy, irhs, 1)
     1               * yf_area_patch ( ixx, iyy, irhs)
     1               / yf_area_patch ( ixx, iyy, ilhs)


      Bxf ( ixx-1, iyy, ilhs, 1) = Bxf ( ixx-1, iyy, irhs, 1)
     1                 * xf_area_patch ( ixx-1, iyy, irhs)
     1                 / xf_area_patch ( ixx-1, iyy, ilhs)

      Byf ( ixx, iyy-1, ilhs, 1) = Byf ( ixx, iyy-1, irhs, 1)
     1                 * yf_area_patch ( ixx, iyy-1, irhs)
     1                 / yf_area_patch ( ixx, iyy-1, ilhs)

! ----------
#else
! ----------

      U ( ixx, iyy, ilhs, :, 1) = U ( ixx, iyy, irhs, :, 1)

      Bxf ( ixx,   iyy, ilhs, 1) = Bxf ( ixx,   iyy, irhs, 1)
      Bxf ( ixx-1, iyy, ilhs, 1) = Bxf ( ixx-1, iyy, irhs, 1)

      Byf ( ixx, iyy,   ilhs, 1) = Byf ( ixx, iyy,   irhs, 1)
      Byf ( ixx, iyy-1, ilhs, 1) = Byf ( ixx, iyy-1, irhs, 1)

      Bzf ( ixx, iyy, ilhs-1, 1) = Bzf ( ixx, iyy, irhs-1, 1)

! ----------
#endif
! ----------

! --------------------

      END IF

! ------------------------------

      END IF

! ------------------------------

      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO izz = 1, NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ilhs, irhs,                  &
!$OMP&    rlen, xu, yu, zu, rho_vr)                                     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, ob_8pi,         &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$OMP&    v_volume_patch, v_centroid_patch,                             &
!$OMP&    imask_active_patch, bcarr_r, U, Bxf, Byf, Bzf)                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$acc&    v_volume_patch, v_centroid_patch,                             &
!$acc&   ob_8pi, imask_active_patch, bcarr_r, U, Bxf, Byf, Bzf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ilhs, irhs)

      DO izz = 1, NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

! ------------------------------

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ------------------------------
! Outflow/Continuitive at upper r-boundary.

      IF ( bcarr_r ( ixx, iyy, 2) == 2 ) THEN

! --------------------

      ilhs = nr_zones + izz
      irhs = nr_zones

      U ( ixx, iyy, ilhs, :, 1) = U ( ixx, iyy, irhs, :, 1)
     1             * v_volume_patch ( ixx, iyy, irhs)
     1             / v_volume_patch ( ixx, iyy, ilhs)

      Bzf ( ixx, iyy, ilhs, 1) = Bzf ( ixx, iyy, irhs, 1)
     1               * zf_area_patch ( ixx, iyy, irhs)
     1               / zf_area_patch ( ixx, iyy, ilhs)

      Bxf ( ixx, iyy, ilhs, 1) = Bxf ( ixx, iyy, irhs, 1)
     1               * xf_area_patch ( ixx, iyy, irhs)
     1               / xf_area_patch ( ixx, iyy, ilhs)

      Byf ( ixx, iyy, ilhs, 1) = Byf ( ixx, iyy, irhs, 1)
     1               * yf_area_patch ( ixx, iyy, irhs)
     1               / yf_area_patch ( ixx, iyy, ilhs)


      Bxf ( ixx-1, iyy, ilhs, 1) = Bxf ( ixx-1, iyy, irhs, 1)
     1                 * xf_area_patch ( ixx-1, iyy, irhs)
     1                 / xf_area_patch ( ixx-1, iyy, ilhs)

      Byf ( ixx, iyy-1, ilhs, 1) = Byf ( ixx, iyy-1, irhs, 1)
     1                 * yf_area_patch ( ixx, iyy-1, irhs)
     1                 / yf_area_patch ( ixx, iyy-1, ilhs)

! --------------------
! Maintains positive radial velocity at outer radial boundary.
#if( CHOOSE_APPL_CAK == 1)

        rlen = SQRT ( v_centroid_patch ( ixx, iyy, irhs, 1)**2
     1              + v_centroid_patch ( ixx, iyy, irhs, 2)**2
     1              + v_centroid_patch ( ixx, iyy, irhs, 3)**2)

        xu = v_centroid_patch ( ixx, iyy, irhs, 1) / rlen
        yu = v_centroid_patch ( ixx, iyy, irhs, 2) / rlen
        zu = v_centroid_patch ( ixx, iyy, irhs, 3) / rlen

        rho_vr = U ( ixx, iyy, irhs, 2, 1) * xu
     1         + U ( ixx, iyy, irhs, 3, 1) * yu
     1         + U ( ixx, iyy, irhs, 4, 1) * zu

! ----------
        IF( rho_vr < 0.0 ) THEN
! ----------

        U ( ixx, iyy, ilhs, 5, 1) = U ( ixx, iyy, ilhs, 5, 1)
     1 - 0.5 * ( U ( ixx, iyy, ilhs, 2, 1)**2
     1         + U ( ixx, iyy, ilhs, 3, 1)**2
     1         + U ( ixx, iyy, ilhs, 4, 1)**2) / U ( ixx, iyy, ilhs,1,1)

        rho_vr = MAX ( 0.0, rho_vr)

        U ( ixx, iyy, ilhs, 2, 1) = rho_vr * xu
        U ( ixx, iyy, ilhs, 3, 1) = rho_vr * yu
        U ( ixx, iyy, ilhs, 4, 1) = rho_vr * zu

        U ( ixx, iyy, ilhs, 5, 1) = U ( ixx, iyy, ilhs, 5, 1)
     1 + 0.5 * ( U ( ixx, iyy, ilhs, 2, 1)**2
     1         + U ( ixx, iyy, ilhs, 3, 1)**2
     1         + U ( ixx, iyy, ilhs, 4, 1)**2) / U ( ixx, iyy, ilhs,1,1)

! ----------
        ENDIF
! ----------

#endif

! --------------------

#if( ( USE_B0B1_FORMULATION == 1 ) && ( SET_0_OUTER_B1 == 1))

        U ( ixx, iyy, ilhs, 5, 1) = U ( ixx, iyy, ilhs, 5, 1)
     1    - ob_8pi * ( U ( ixx, iyy, ilhs, 6, 1)**2
     1               + U ( ixx, iyy, ilhs, 7, 1)**2
     1               + U ( ixx, iyy, ilhs, 8, 1)**2)

        U ( ixx, iyy, ilhs, 6, 1) = 0.0
        U ( ixx, iyy, ilhs, 7, 1) = 0.0
        U ( ixx, iyy, ilhs, 8, 1) = 0.0

        Bzf ( ixx, iyy, ilhs, 1) = 0.0

! setting the Bzf for the last zone to zero. Note: irhs == nr_zones.
! Commenting out because it makes the divB error non-zero!
!     Bzf ( ixx, iyy, irhs, 1) = 0.0

      Bxf ( ixx,   iyy, ilhs, 1) = 0.0
      Bxf ( ixx-1, iyy, ilhs, 1) = 0.0

      Byf ( ixx, iyy,   ilhs, 1) = 0.0
      Byf ( ixx, iyy-1, ilhs, 1) = 0.0

#endif

! ------------------------------

      END IF  ! End of "IF ( bcarr_r ( ixx, iyy, 2) == 2 ) THEN".

! ------------------------------

      END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

! ------------------------------

      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO izz = 1, NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

      END SUBROUTINE PAD_BOUNDARY_R
