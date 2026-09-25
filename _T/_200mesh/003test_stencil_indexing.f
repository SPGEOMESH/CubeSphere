




      PROGRAM TEST_STENCIL_INDEXING

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include <cubesphere/000com/riemann.com>

      INTEGER ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,
     1        ix_patch, iy_patch, ixx, iyy,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      PARAMETER ( ntp_zones = 8, nr_zones = 2, NGhost = 4)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------

      OPEN ( UNIT = 6, FILE = "test_stencil.out")

      ntp_patches = 3
      nr_patches = 1

      ix_patch = 2
      iy_patch = 2

! "imask_active_patch, istencil_patch" are copied to GPU.
! All scalars guaranteed to be copied.

!$acc enter data copyin( imask_active_patch, istencil_patch )

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

!$acc update host( imask_active_patch, istencil_patch )

! "imask_active_patch, istencil_patch" are brought back to CPU from GPU.

! ----------------------------------------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1   "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU = "
      WRITE ( 6, *)
     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------
! Print, Down - along x-axis, iyy = 1.

      iyy = 1

      WRITE ( 6, *)"Down - along x-axis, iyy = 1"
      WRITE ( 6, *)" "

      DO ixx = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 4)",
     1              istencil_patch( ixx, iyy, 1, 4),
     1              istencil_patch( ixx, iyy, 2, 4)

      END DO

! ----------------------------------------------------------------------
! Print, Up - along x-axis, iyy = ntp_zones.

      iyy = ntp_zones

      WRITE ( 6, *)"Up - along x-axis, iyy = ntp_zones"
      WRITE ( 6, *)" "

      DO ixx = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 3) = ",
     1              istencil_patch( ixx, iyy, 1, 3),
     1              istencil_patch( ixx, iyy, 2, 3)

      END DO


! ----------------------------------------------------------------------
! Print, Left - along y-axis, ixx = 1.

      ixx = 1

      WRITE ( 6, *)"Left - along y-axis, ixx = 1"
      WRITE ( 6, *)" "

      DO iyy = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 2)",
     1              istencil_patch( ixx, iyy, 1, 2),
     1              istencil_patch( ixx, iyy, 2, 2)

      END DO

! ----------------------------------------------------------------------
! Print, Right - along y-axis, ixx = ntp_zones.

      ixx = ntp_zones

      WRITE ( 6, *)"Right - along y-axis, ixx = ntp_zones"
      WRITE ( 6, *)" "

      DO iyy = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 1)",
     1              istencil_patch( ixx, iyy, 1, 1),
     1              istencil_patch( ixx, iyy, 2, 1)

      END DO

! ----------------------------------------------------------------------

      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      END PROGRAM TEST_STENCIL_INDEXING
