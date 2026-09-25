




      SUBROUTINE STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! If most zones "ixx + 1" and "ixx - 1" gives the zones to the Right or Left.
! Similarly, "iyy + 1" and "iyy - 1" gives the Up or Down zones.
! However, at the tri-corners of a cube sphere mesh, this ceases to be the
! case. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

! The stencils are wide enough to permit construction of second and third order
! schemes. They can even be extended to WENO-AO for fourth order.

! The zone nomenclature, arrangement and indexing is as follows:-

!                   UU                                   7

!              LU   U   RU                         10    3     9

!        LL    L    C    R    RR     -->     6     2     0     1     5

!              LD   D   RD                         11    4    12

!                   DD                                   8

! "ntp_zones, nr_zones, NGhost" are the number of zones in the theta-phi direction,
! the radial direction and the ghost zones on each patch. Any patch in the entire
! AMR hierarchy has this fixed number of zones and ghost zones.
! This is baked into the data structures for the whole code.

! "ntp_patches, nr_patches" are the total number of patches in the theta-phi and
! radial directions IF that level were indeed completely full. Other than the
! base level, the AMR levels do not have to be completely full.

! "ix_patch, iy_patch" is the patch at level "ilevel" that we indeed want.
! Again, this is based on expecting a full level whether there are physical levels
! there or not.

! "imask_active_patch" (input) is 1 for active zones, 0 for inactive ones
! that may lie close to the corners.

! "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU" (output) labels the
! four corners in the xy-plane. 1 if this is a tri-corner patch; 0 if this
! is a normal quad-corner patch.

! "istencil_patch ( :, :, 2, 12)" (output)
! "istencil_patch ( :, :, 2, 1:4)" gives us Right, Left, Up, Down zones in xy-plane.
! "istencil_patch ( :, :, 2, 5:8)" gives us the zones that are two ones over in the
! Right, Left, Up, Down directions.
! "istencil_patch ( :, :, 2, 9:12)" gives us the RU, LU, LD, RD zones.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,
     1        ix_patch, iy_patch,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

! Output variable, mask which is zone-centered and 1 if zone is active, 0 if inactive.
      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------
! Variables below this line are internal to the subroutine.

      INTEGER ixx, iyy, ixxp1, iyyp1, ixxm1, iyym1,
     2        ixxp2, iyyp2, ixxm2, iyym2

! ----------------------------------------------------------------------

! By default, set "itri_patch_LD:4" to 0 -- i.e. it is a normal quad-corner patch.
! Then identify tri-corner patches based on analyzing each of the four corners
! of the patch in the xy-plane.

      itri_patch_LD = 0  ;  itri_patch_RD = 0  ! Default to all quad-corner patches.
      itri_patch_RU = 0  ;  itri_patch_LU = 0


! Left-Down vertex of this patch is a tri-corner.
      IF ( ( ix_patch .EQ. 1) .AND.
     1     ( iy_patch .EQ. 1) ) itri_patch_LD = 1

! Right-Down vertex of this patch is a tri-corner.
      IF ( ( ix_patch .EQ. ntp_patches) .AND.
     1     ( iy_patch .EQ. 1) ) itri_patch_RD = 1

! Right-Up vertex of this patch is a tri-corner.
      IF ( ( ix_patch .EQ. ntp_patches) .AND.
     1     ( iy_patch .EQ. ntp_patches) ) itri_patch_RU = 1

! Left-Up vertex of this patch is a tri-corner.
      IF ( ( ix_patch .EQ. 1) .AND.
     1     ( iy_patch .EQ. ntp_patches) ) itri_patch_LU = 1



!$acc update device ( itri_patch_LD, itri_patch_RD, itri_patch_RU,      &
!$acc&                itri_patch_LU)

! ----------------------------------------------------------------------

! Now with the assumption that all corners are quad-corners, fill in
! "istencil_patch ( :, :, 2, 12)"

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, ixxp1, iyyp1, ixxp2, iyyp2,       &
!$OMP&    ixxm1, iyym1, ixxm2, iyym2)                                   &
!$OMP&            SHARED ( ntp_zones, NGhost, istencil_patch)           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel if ( do_on_gpu) vector_length(64)                        &
!$acc& present( istencil_patch)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, ixxp1, iyyp1, ixxp2, iyyp2,               &
!$acc&    ixxm1, iyym1, ixxm2, iyym2)

      DO ixx = 1 - NGhost, ntp_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------

        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)

        ixxp2 = MIN0 ( ixx + 2, ntp_zones + NGhost)
        iyyp2 = MIN0 ( iyy + 2, ntp_zones + NGhost)

        ixxm1 = MAX0 ( ixx - 1, 1 - NGhost)
        iyym1 = MAX0 ( iyy - 1, 1 - NGhost)

        ixxm2 = MAX0 ( ixx - 2, 1 - NGhost)
        iyym2 = MAX0 ( iyy - 2, 1 - NGhost)

! ----------

        istencil_patch ( ixx, iyy, 1, 1) = ixxp1  ! Zone that is one to the Right.
        istencil_patch ( ixx, iyy, 2, 1) = iyy

        istencil_patch ( ixx, iyy, 1, 2) = ixxm1  ! Zone that is one to the Left.
        istencil_patch ( ixx, iyy, 2, 2) = iyy

        istencil_patch ( ixx, iyy, 1, 3) = ixx  ! Zone that is one to the Up.
        istencil_patch ( ixx, iyy, 2, 3) = iyyp1

        istencil_patch ( ixx, iyy, 1, 4) = ixx  ! Zone that is one to the Down.
        istencil_patch ( ixx, iyy, 2, 4) = iyym1

! ----------

        istencil_patch ( ixx, iyy, 1, 5) = ixxp2  ! Zone that is two to the Right.
        istencil_patch ( ixx, iyy, 2, 5) = iyy

        istencil_patch ( ixx, iyy, 1, 6) = ixxm2  ! Zone that is two to the Left.
        istencil_patch ( ixx, iyy, 2, 6) = iyy

        istencil_patch ( ixx, iyy, 1, 7) = ixx  ! Zone that is two to the Up.
        istencil_patch ( ixx, iyy, 2, 7) = iyyp2

        istencil_patch ( ixx, iyy, 1, 8) = ixx  ! Zone that is two to the Down.
        istencil_patch ( ixx, iyy, 2, 8) = iyym2

! ----------

        istencil_patch ( ixx, iyy, 1, 9) = ixxp1  ! Zone that is Right-Up.
        istencil_patch ( ixx, iyy, 2, 9) = iyyp1

        istencil_patch ( ixx, iyy, 1, 10) = ixxm1  ! Zone that is Left-Up.
        istencil_patch ( ixx, iyy, 2, 10) = iyyp1

        istencil_patch ( ixx, iyy, 1, 11) = ixxm1  ! Zone that is Left-Down.
        istencil_patch ( ixx, iyy, 2, 11) = iyym1

        istencil_patch ( ixx, iyy, 1, 12) = ixxp1  ! Zone that is Right-Down.
        istencil_patch ( ixx, iyy, 2, 12) = iyym1

! ----------

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Left-Down vertex is missing one quadrant.

! The "( 0, 1)" zone will be missing a RD zone.
! The "( 1, 0)" zone will be missing a LU zone.

!$acc kernels if ( do_on_gpu)

      IF ( itri_patch_LD .EQ. 1) THEN

! ------------------------------

      iyy = 1

!$acc loop vector
      DO ixx = 1 - NGhost, 0

! ----------

! Zone that is one to the Down.
        istencil_patch ( ixx, iyy, 1, 4) = 1
        istencil_patch ( ixx, iyy, 2, 4) = ixx

! Zone that is two to the Down.
        istencil_patch ( ixx, iyy, 1, 8) = 2
        istencil_patch ( ixx, iyy, 2, 8) = ixx

! ----------

! Zone that is Left-Down.
        istencil_patch ( ixx, iyy, 1, 11) = 1
        istencil_patch ( ixx, iyy, 2, 11) = MAX0 ( ixx - 1, 1 - NGhost)

! Zone that is Right-Down.
        istencil_patch ( ixx, iyy, 1, 12) = 1
        istencil_patch ( ixx, iyy, 2, 12) = ixx + 1

        IF ( ixx .EQ. 0) THEN
! Zone that is two Down and one Left instead, because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 12) = 2
        istencil_patch ( ixx, iyy, 2, 12) = ixx - 1
        END IF

! ----------

! Zone that is two to the Down; one layer away.
        istencil_patch ( ixx, iyy + 1, 1, 8) = 1
        istencil_patch ( ixx, iyy + 1, 2, 8) = ixx

! ----------

      END DO

! ------------------------------

      ixx = 1

!$acc loop vector
      DO iyy = 1 - NGhost, 0

! ----------

! Zone that is one to the Left.
        istencil_patch ( ixx, iyy, 1, 2) = iyy
        istencil_patch ( ixx, iyy, 2, 2) = 1

! Zone that is two to the Left.
        istencil_patch ( ixx, iyy, 1, 6) = iyy
        istencil_patch ( ixx, iyy, 2, 6) = 2

! ----------

! Zone that is Left-Up.
        istencil_patch ( ixx, iyy, 1, 10) = iyy + 1
        istencil_patch ( ixx, iyy, 2, 10) = 1

        IF ( iyy .EQ. 0) THEN
! Zone that is two Left and one Down instead, because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 10) = iyy - 1
        istencil_patch ( ixx, iyy, 2, 10) = 2
        END IF

! Zone that is Left-Down.
        istencil_patch ( ixx, iyy, 1, 11) = MAX0 ( iyy - 1, 1 - NGhost)
        istencil_patch ( ixx, iyy, 2, 11) = 1

! ----------

! Zone that is two to the Left; one layer away.
        istencil_patch ( ixx + 1, iyy, 1, 6) = iyy
        istencil_patch ( ixx + 1, iyy, 2, 6) = 1

! ----------

      END DO

! ------------------------------

! Zone that is two Down and one Right of the vertex, because it too will give
! a unique stencil.

        istencil_patch ( 1, 1, 1, 11) = 2
        istencil_patch ( 1, 1, 2, 11) = -1

! ------------------------------

      END IF  ! End of "IF ( itri_patch_LD .EQ. 1) THEN".

! ----------------------------------------------------------------------
! Right-Down vertex is missing one quadrant.

! The "( ntp_zones + 1, 1)" zone will be missing a LD zone.
! The "( ntp_zones, 0)" zone will be missing a RU zone.

      IF ( itri_patch_RD .EQ. 1) THEN

! ------------------------------

      iyy = 1

!$acc loop vector
      DO ixx = ntp_zones + 1, ntp_zones + NGhost

! ----------

! Zone that is one to the Down.
        istencil_patch ( ixx, iyy, 1, 4) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 4) = ntp_zones - ixx + 1

! Zone that is two to the Down.
        istencil_patch ( ixx, iyy, 1, 8) = ntp_zones - 1
        istencil_patch ( ixx, iyy, 2, 8) = ntp_zones - ixx + 1

! ----------

! Zone that is Left-Down.
        istencil_patch ( ixx, iyy, 1, 11) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 11) = ntp_zones - ixx + 2

        IF ( ixx .EQ. ntp_zones + 1) THEN
! Zone that is two Down and one Right, because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 11) = ntp_zones - 1
        istencil_patch ( ixx, iyy, 2, 11) = ntp_zones - ixx
        END IF

! Zone that is Right-Down.
        istencil_patch ( ixx, iyy, 1, 12) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 12) =
     1                             MAX0 ( ntp_zones - ixx, 1 - NGhost)

! ----------

! Zone that is two to the Down; one layer away.
        istencil_patch ( ixx, iyy + 1, 1, 8) = ntp_zones
        istencil_patch ( ixx, iyy + 1, 2, 8) = ntp_zones - ixx + 1

! ----------

      END DO

! ------------------------------

      ixx = ntp_zones

!$acc loop vector
      DO iyy = 1 - NGhost, 0

! ----------

! Zone that is one to the Right.
        istencil_patch ( ixx, iyy, 1, 1) = ntp_zones + 1 - iyy
        istencil_patch ( ixx, iyy, 2, 1) = 1

! Zone that is two to the Right.
        istencil_patch ( ixx, iyy, 1, 5) = ntp_zones + 1 - iyy
        istencil_patch ( ixx, iyy, 2, 5) = 2

! ----------

! Zone that is Right-Up.
        istencil_patch ( ixx, iyy, 1, 9) = ntp_zones - iyy
        istencil_patch ( ixx, iyy, 2, 9) = 1

        IF ( iyy .EQ. 0) THEN
! Zone that is two to the Right and one Down, because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 9) = ntp_zones + 2 - iyy
        istencil_patch ( ixx, iyy, 2, 9) = 2
        END IF

! Zone that is Right-Down.
        istencil_patch ( ixx, iyy, 1, 12) =
     1                   MIN0 ( ntp_zones + 2 - iyy, ntp_zones + NGhost)
        istencil_patch ( ixx, iyy, 2, 12) = 1

! ----------

! Zone that is two to the Right; one layer away.
        istencil_patch ( ixx - 1, iyy, 1, 5) = ntp_zones + 1 - iyy
        istencil_patch ( ixx - 1, iyy, 2, 5) = 1

! ----------

      END DO

! ------------------------------

! Zone that is two Down and one Left of the vertex, because it too will give
! a unique stencil.

        istencil_patch ( ntp_zones, 1, 1, 12) = ntp_zones - 1
        istencil_patch ( ntp_zones, 1, 2, 12) = -1

! ------------------------------

      END IF  ! End of "IF ( itri_patch_RD .EQ. 1) THEN".

! ----------------------------------------------------------------------
! Right-Up vertex is missing one quadrant.

! The "( ntp_zones + 1, ntp_zones)" zone will be missing a LU zone.
! The "( ntp_zones, ntp_zones + 1)" zone will be missing a RD zone.

      IF ( itri_patch_RU .EQ. 1) THEN

! ------------------------------

      iyy = ntp_zones

!$acc loop vector
      DO ixx = ntp_zones + 1, ntp_zones + NGhost

! ----------

! Zone that is one to the Up.
        istencil_patch ( ixx, iyy, 1, 3) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 3) = ixx

! Zone that is two to the Up.
        istencil_patch ( ixx, iyy, 1, 7) = ntp_zones - 1
        istencil_patch ( ixx, iyy, 2, 7) = ixx

! ----------

! Zone that is Left-Up.
        istencil_patch ( ixx, iyy, 1, 10) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 10) = ixx - 1

        IF ( ixx .EQ. ntp_zones + 1) THEN
! Zone that is two to the Up and one Right; because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 10) = ntp_zones - 1
        istencil_patch ( ixx, iyy, 2, 10) = ixx + 1
        END IF

! Zone that is Right-Up.
        istencil_patch ( ixx, iyy, 1, 9) = ntp_zones
        istencil_patch ( ixx, iyy, 2, 9) =
     1                   MIN0 ( ixx + 1, ntp_zones + NGhost)

! ----------

! Zone that is two to the Up; one layer away.
        istencil_patch ( ixx, iyy - 1, 1, 7) = ntp_zones
        istencil_patch ( ixx, iyy - 1, 2, 7) = ixx

! ----------

      END DO

! ------------------------------

      ixx = ntp_zones

!$acc loop vector
      DO iyy = ntp_zones + 1, ntp_zones + NGhost

! ----------

! Zone that is one to the Right.
        istencil_patch ( ixx, iyy, 1, 1) = iyy
        istencil_patch ( ixx, iyy, 2, 1) = ntp_zones

! Zone that is two to the Right.
        istencil_patch ( ixx, iyy, 1, 5) = iyy
        istencil_patch ( ixx, iyy, 2, 5) = ntp_zones - 1

! ----------

! Zone that is Right-Down.
        istencil_patch ( ixx, iyy, 1, 12) = iyy - 1
        istencil_patch ( ixx, iyy, 2, 12) = ntp_zones

        IF ( iyy .EQ. ntp_zones + 1) THEN
! Zone that is two Right and one Up; because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 12) = iyy + 1
        istencil_patch ( ixx, iyy, 2, 12) = ntp_zones - 1
        END IF

! Zone that is Right-Up.
        istencil_patch ( ixx, iyy, 1, 9) =
     1                   MIN0 ( iyy + 1, ntp_zones + NGhost)
        istencil_patch ( ixx, iyy, 2, 9) = ntp_zones

! ----------

! Zone that is two to the Right; one layer away.
        istencil_patch ( ixx - 1, iyy, 1, 5) = iyy
        istencil_patch ( ixx - 1, iyy, 2, 5) = ntp_zones

! ----------

      END DO

! ------------------------------

! Zone that is two Up and one Left of the vertex; because it too will give
! a unique stencil.

        istencil_patch ( ntp_zones, ntp_zones, 1, 9) = ntp_zones - 1
        istencil_patch ( ntp_zones, ntp_zones, 2, 9) = ntp_zones + 2

! ------------------------------

      END IF  ! End of "IF ( itri_patch_RU .EQ. 1) THEN".

! ----------------------------------------------------------------------
! Left-Up vertex is missing one quadrant.

! The "( 1, ntp_zones + 1)" zone will be missing a LD zone.
! The "( 0, ntp_zones)" zone will be missing a RU zone.

      IF ( itri_patch_LU .EQ. 1) THEN

! ------------------------------

      ixx = 1

!$acc loop vector
      DO iyy = ntp_zones + 1, ntp_zones + NGhost

! ----------

! Zone that is one to the Left.
        istencil_patch ( ixx, iyy, 1, 2) = ntp_zones - iyy + 1
        istencil_patch ( ixx, iyy, 2, 2) = ntp_zones

! Zone that is two to the Left.
        istencil_patch ( ixx, iyy, 1, 6) = ntp_zones - iyy + 1
        istencil_patch ( ixx, iyy, 2, 6) = ntp_zones - 1

! ----------

! Zone that is Left-Down.
        istencil_patch ( ixx, iyy, 1, 11) = ntp_zones - iyy + 2
        istencil_patch ( ixx, iyy, 2, 11) = ntp_zones

        IF ( iyy .EQ. ntp_zones + 1) THEN
! Zone that is two Left and one Up; because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 11) = ntp_zones - iyy
        istencil_patch ( ixx, iyy, 2, 11) = ntp_zones - 1
        END IF

! Zone that is Left-Up.
        istencil_patch ( ixx, iyy, 1, 10) =
     1                   MAX0 ( ntp_zones - iyy, 1 - NGhost)
        istencil_patch ( ixx, iyy, 2, 10) = ntp_zones

! ----------

! Zone that is two to the Left; one layer away.
        istencil_patch ( ixx + 1, iyy, 1, 6) = ntp_zones - iyy + 1
        istencil_patch ( ixx + 1, iyy, 2, 6) = ntp_zones

! ----------

      END DO

! ------------------------------

      iyy = ntp_zones

!$acc loop vector
      DO ixx = 1 - NGhost, 0

! ----------

! Zone that is one to the Up.
        istencil_patch ( ixx, iyy, 1, 3) = 1
        istencil_patch ( ixx, iyy, 2, 3) = ntp_zones - ixx + 1

! Zone that is two to the Up.
        istencil_patch ( ixx, iyy, 1, 7) = 2
        istencil_patch ( ixx, iyy, 2, 7) = ntp_zones - ixx + 1

! ----------

! Zone that is Right-Up.
        istencil_patch ( ixx, iyy, 1, 9) = 1
        istencil_patch ( ixx, iyy, 2, 9) = ntp_zones - ixx

        IF ( ixx .EQ. 0) THEN
! Zone that is two Up and one Left; because it too will give
! a unique stencil.
        istencil_patch ( ixx, iyy, 1, 9) = 2
        istencil_patch ( ixx, iyy, 2, 9) = ntp_zones - ixx + 2
        END IF

! Zone that is Left-Up.
        istencil_patch ( ixx, iyy, 1, 10) = 1
        istencil_patch ( ixx, iyy, 2, 10) =
     1                   MIN0 ( ntp_zones - ixx + 2, ntp_zones + NGhost)

! ----------

! Zone that is two to the Up; one layer away.
        istencil_patch ( ixx, iyy - 1, 1, 7) = 1
        istencil_patch ( ixx, iyy - 1, 2, 7) = ntp_zones - ixx + 1

! ----------

      END DO

! ------------------------------

! Zone that is two Up and one Right of the vertex; because it too will give
! a unique stencil.

        istencil_patch ( 1, ntp_zones, 1, 10) = 2
        istencil_patch ( 1, ntp_zones, 2, 10) = ntp_zones + 2

! ------------------------------

      END IF  ! End of "IF ( itri_patch_LU .EQ. 1) THEN".

!$acc end kernels

! ----------------------------------------------------------------------

      END SUBROUTINE STENCIL_INDEXING_TRICORNERS




