





      SUBROUTINE COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara and Sethu Subramanyan

! This subroutine takes in "vert_coords_patch" which contains the coordinates
! of the vertices of each patch of the cube sphere mesh. It also takes in
! "imask_active_patch" which is 1 for good zones, 0 for bad ones. It then performs
! isoparametric mapping of degree two in the volume of the mesh.

! It returns the locations of the quadrature points in the volume :-
! "v_quadpts_patch".

! It also returns the weights associated with the quadrature points in the volume :-
! "v_quadwt_patch".

! It also returns the derivative of the position vector relative to reference
! coordinates, i.e. the tangent vectors at the quadrature points in the volume :-
! "v_quadtangent_patch".

! It also returns the volume of each zone :- "v_volume_patch".

! It also returns the centroid of each zone :- "v_centroid_patch".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost

! Input variable, gives us the vertices of a cube sphere mesh.
      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

! Input variable, mask which is zone-centered and 1 if zone is active, 0 if inactive.
      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the scalar weight  for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the derivative of the position vector relative to
! reference coordinates for each of the "NUMORDER * NUMORDER * NUMORDER" number of
! quadrature points in the volume. These are needed in the ADER update step.
      REAL v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the volume of the zone.
      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! Output variable, gives us the centroid of the zone.
      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------------------------------------------------------------------
! Variables internal to this subroutine defined here.

      INTEGER ixx, iyy, izz, ixxm1, iyym1, izzm1, i

      REAL, DIMENSION ( NUMORDER) :: xg, yg, zg, lin_wt

      REAL, DIMENSION ( 20, NUMORDER * NUMORDER * NUMORDER) :: basis_3D,
     1      d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz

      REAL x, y, z, q11 ( 1), vec_in ( 20)

      REAL, DIMENSION ( NUMORDER * NUMORDER * NUMORDER) :: wt, vec_out
      REAL, DIMENSION ( 1, NUMORDER * NUMORDER * NUMORDER) :: qwt

! ----------

      REAL vert ( 20, 3)

      REAL radius1, radius2, tempaa, tempab

      REAL, DIMENSION ( NUMORDER * NUMORDER * NUMORDER) ::
     1      Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, Det

      REAL v_centroid_temp ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)
!$acc data create( v_centroid_temp )

! ----------------------------------------------------------------------
! STEP 1.a) Build coordinates and weights of quadrature points in the reference
! element. Gaussian quadrature of successive orders used.

! ------------------------------
#if ( NUMORDER == 2)

      xg ( 1) = - 0.5 * SQRT ( 1.0 / 3.0)
      xg ( 2) = + 0.5 * SQRT ( 1.0 / 3.0)

      yg ( :) = xg ( :)
      zg ( :) = xg ( :)

      lin_wt ( 1) = 0.5
      lin_wt ( 2) = 0.5

#endif
! ------------------------------
#if ( NUMORDER == 3)

      xg(1) = - 0.5 * SQRT ( 3.0 / 5.0)
      xg(2) = 0.0
      xg(3) = 0.5 * SQRT ( 3.0 / 5.0)

      yg ( :) = xg ( :)
      zg ( :) = xg ( :)

      lin_wt ( 1) = 5.0 / 18.0
      lin_wt ( 2) = 8.0 / 18.0
      lin_wt ( 3) = 5.0 / 18.0

#endif
! ------------------------------
#if ( NUMORDER == 4)

      xg(1) = - SQRT( 3.0/7.0 + (2.0/7.0) * SQRT( 6.0/5.0 ) ) / 2.0
      xg(2) = - SQRT( 3.0/7.0 - (2.0/7.0) * SQRT( 6.0/5.0 ) ) / 2.0
      xg(3) = + SQRT( 3.0/7.0 - (2.0/7.0) * SQRT( 6.0/5.0 ) ) / 2.0
      xg(4) = + SQRT( 3.0/7.0 + (2.0/7.0) * SQRT( 6.0/5.0 ) ) / 2.0

      yg ( :) = xg ( :)
      zg ( :) = xg ( :)

      lin_wt ( 1) = ( 18.0 - SQRT( 30.0)) / 72.0
      lin_wt ( 2) = ( 18.0 + SQRT( 30.0)) / 72.0
      lin_wt ( 3) = ( 18.0 + SQRT( 30.0)) / 72.0
      lin_wt ( 4) = ( 18.0 - SQRT( 30.0)) / 72.0

#endif
!-------------------------------
! Weights are symmetric in all 3-directions, so we are free to choose
! the x, y, z loop order, to get the tensor product "wt ( :)".

      i = 0

      DO ixx = 1, NUMORDER
        DO iyy = 1, NUMORDER
          DO izz = 1, NUMORDER

            i = i + 1

            wt ( i) = lin_wt ( ixx) *  lin_wt ( iyy) *  lin_wt ( izz)

          END DO
        END DO
      END DO

! --------------------------------------------------
! STEP 1.b) Now use the above to build the basis functions and their derivatives
! in the 3D reference element. We will use isoparametric mapping with second degree
! polynomial nodal bases. This works at all reconstruction orders.

      i = 0

      DO izz = 1, NUMORDER

        DO iyy = 1, NUMORDER

          DO ixx = 1, NUMORDER

          i = i + 1

          x = xg ( ixx)
          y = yg ( iyy)
          z = zg ( izz)

! --------------------
! Serendipity Basis functions for a 20 Node Cube.

! ----------
! Vertex Nodes - lower z-plane

      basis_3D (  1, i) = (0.5 - x) * (0.5 - y) * (0.5 - z)
     1                              * (-2.0*x - 2.0*y - 2.0*z - 2.0)

      basis_3D (  2, i) = (0.5 + x) * (0.5 - y) * (0.5 - z)
     1                              * ( 2.0*x - 2.0*y - 2.0*z - 2.0)

      basis_3D (  3, i) = (0.5 + x) * (0.5 + y) * (0.5 - z)
     1                              * ( 2.0*x + 2.0*y - 2.0*z - 2.0)

      basis_3D (  4, i) = (0.5 - x) * (0.5 + y) * (0.5 - z)
     1                              * (-2.0*x + 2.0*y - 2.0*z - 2.0)

! ----------
! Vertex Nodes - upper z-plane
      basis_3D (  5, i) = (0.5 - x) * (0.5 - y) * (0.5 + z)
     1                              * (-2.0*x - 2.0*y + 2.0*z - 2.0)

      basis_3D (  6, i) = (0.5 + x) * (0.5 - y) * (0.5 + z)
     1                              * ( 2.0*x - 2.0*y + 2.0*z - 2.0)

      basis_3D (  7, i) = (0.5 + x) * (0.5 + y) * (0.5 + z)
     1                              * ( 2.0*x + 2.0*y + 2.0*z - 2.0)

      basis_3D (  8, i) = (0.5 - x) * (0.5 + y) * (0.5 + z)
     1                              * (-2.0*x + 2.0*y + 2.0*z - 2.0)

! ----------
! x = 0 Nodes
      basis_3D (  9, i) = (1.0 - 4.0*x**2) * (0.5 - y) * (0.5 - z)
      basis_3D ( 10, i) = (1.0 - 4.0*x**2) * (0.5 + y) * (0.5 - z)
      basis_3D ( 11, i) = (1.0 - 4.0*x**2) * (0.5 + y) * (0.5 + z)
      basis_3D ( 12, i) = (1.0 - 4.0*x**2) * (0.5 - y) * (0.5 + z)

! ----------
! y = 0 Nodes
      basis_3D ( 13, i) = (0.5 - x) * (1.0 - 4.0*y**2) * (0.5 - z)
      basis_3D ( 14, i) = (0.5 - x) * (1.0 - 4.0*y**2) * (0.5 + z)
      basis_3D ( 15, i) = (0.5 + x) * (1.0 - 4.0*y**2) * (0.5 + z)
      basis_3D ( 16, i) = (0.5 + x) * (1.0 - 4.0*y**2) * (0.5 - z)

! ----------
! z = 0 Nodes
      basis_3D ( 17, i) = (0.5 - x) * (0.5 - y) * (1.0 - 4.0*z**2)
      basis_3D ( 18, i) = (0.5 + x) * (0.5 - y) * (1.0 - 4.0*z**2)
      basis_3D ( 19, i) = (0.5 + x) * (0.5 + y) * (1.0 - 4.0*z**2)
      basis_3D ( 20, i) = (0.5 - x) * (0.5 + y) * (1.0 - 4.0*z**2)

! --------------------
! x - derivative of Serendipity Basis functions for a 20 Node Cube.

! ----------
! Vertex Nodes - lower z-plane
      d_basis_3D_dx (  1, i) =  ( (-1.0 + 2.0*y) * (-1.0 + 2.0*z)
     1                  * ( 1.0 + 4.0*x + 2.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dx (  2, i) = -( (-1.0 + 2.0*y) * (-1.0 + 2.0*z)
     1                  * ( 1.0 - 4.0*x + 2.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dx (  3, i) = -( ( 1.0 + 2.0*y) * (-1.0 + 2.0*z)
     1                  * (-1.0 + 4.0*x + 2.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dx (  4, i) = -( ( 1.0 + 2.0*y) * (-1.0 + 2.0*z)
     1                  * ( 1.0 + 4.0*x - 2.0*y + 2.0*z) ) * 0.25

! ----------
! Vertex Nodes - upper z-plane
      d_basis_3D_dx (  5, i) = -( (-1.0 + 2.0*y) * ( 1.0 + 2.0*z)
     1                  * ( 1.0 + 4.0*x + 2.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dx (  6, i) =  ( (-1.0 + 2.0*y) * ( 1.0 + 2.0*z)
     1                  * ( 1.0 - 4.0*x + 2.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dx (  7, i) =  ( ( 1.0 + 2.0*y) * ( 1.0 + 2.0*z)
     1                  * (-1.0 + 4.0*x + 2.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dx (  8, i) =  ( ( 1.0 + 2.0*y) * ( 1.0 + 2.0*z)
     1                  * ( 1.0 + 4.0*x - 2.0*y - 2.0*z) ) * 0.25

! ----------
! x = 0 Nodes
      d_basis_3D_dx (  9, i) = -2.0*x * (-1.0 + 2.0*y) * (-1.0 + 2.0*z)
      d_basis_3D_dx ( 10, i) =  2.0*x * ( 1.0 + 2.0*y) * (-1.0 + 2.0*z)
      d_basis_3D_dx ( 11, i) = -2.0*x * ( 1.0 + 2.0*y) * ( 1.0 + 2.0*z)
      d_basis_3D_dx ( 12, i) =  2.0*x * (-1.0 + 2.0*y) * ( 1.0 + 2.0*z)

! ----------
! y = 0 Nodes
      d_basis_3D_dx ( 13, i) = -( (1.0 - 4.0*y**2) * (0.5 - z) )
      d_basis_3D_dx ( 14, i) = -( (1.0 - 4.0*y**2) * (0.5 + z) )
      d_basis_3D_dx ( 15, i) =  ( (1.0 - 4.0*y**2) * (0.5 + z) )
      d_basis_3D_dx ( 16, i) =  ( (1.0 - 4.0*y**2) * (0.5 - z) )

! ----------
! z = 0 Nodes
      d_basis_3D_dx ( 17, i) = -( (0.5 - y) * (1.0 - 4.0*z**2) )
      d_basis_3D_dx ( 18, i) =  ( (0.5 - y) * (1.0 - 4.0*z**2) )
      d_basis_3D_dx ( 19, i) =  ( (0.5 + y) * (1.0 - 4.0*z**2) )
      d_basis_3D_dx ( 20, i) = -( (0.5 + y) * (1.0 - 4.0*z**2) )


! --------------------
! y - derivative of Serendipity Basis functions for a 20 Node Cube.

! ----------
! Vertex Nodes - lower z-plane
      d_basis_3D_dy (  1, i) =  ( (-1.0 + 2.0*x) * (-1.0 + 2.0*z)
     1                  * ( 1.0 + 2.0*x + 4.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dy (  2, i) =  ( ( 1.0 + 2.0*x) * (-1.0 + 2.0*z)
     1                  * (-1.0 + 2.0*x - 4.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dy (  3, i) = -( ( 1.0 + 2.0*x) * (-1.0 + 2.0*z)
     1                  * (-1.0 + 2.0*x + 4.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dy (  4, i) = -( (-1.0 + 2.0*x) * (-1.0 + 2.0*z)
     1                  * ( 1.0 + 2.0*x - 4.0*y + 2.0*z) ) * 0.25

! ----------
! Vertex Nodes - upper z-plane
      d_basis_3D_dy (  5, i) = -( (-1.0 + 2.0*x) * ( 1.0 + 2.0*z)
     1                  * ( 1.0 + 2.0*x + 4.0*y - 2.0*z) ) * 0.25

      d_basis_3D_dy (  6, i) = -( ( 1.0 + 2.0*x) * ( 1.0 + 2.0*z)
     1                  * (-1.0 + 2.0*x - 4.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dy (  7, i) =  ( ( 1.0 + 2.0*x) * ( 1.0 + 2.0*z)
     1                  * (-1.0 + 2.0*x + 4.0*y + 2.0*z) ) * 0.25

      d_basis_3D_dy (  8, i) =  ( (-1.0 + 2.0*x) * ( 1.0 + 2.0*z)
     1                  * ( 1.0 + 2.0*x - 4.0*y - 2.0*z) ) * 0.25

! ----------
! x = 0 Nodes
      d_basis_3D_dy (  9, i) = -( (1.0 - 4.0*x**2) * (0.5 - z) )
      d_basis_3D_dy ( 10, i) =  ( (1.0 - 4.0*x**2) * (0.5 - z) )
      d_basis_3D_dy ( 11, i) =  ( (1.0 - 4.0*x**2) * (0.5 + z) )
      d_basis_3D_dy ( 12, i) = -( (1.0 - 4.0*x**2) * (0.5 + z) )

! ----------
! y = 0 Nodes
      d_basis_3D_dy ( 13, i) = -8.0 * (0.5 - x) * y * (0.5 - z)
      d_basis_3D_dy ( 14, i) = -8.0 * (0.5 - x) * y * (0.5 + z)
      d_basis_3D_dy ( 15, i) = -8.0 * (0.5 + x) * y * (0.5 + z)
      d_basis_3D_dy ( 16, i) = -8.0 * (0.5 + x) * y * (0.5 - z)

! ----------
! z = 0 Nodes
      d_basis_3D_dy ( 17, i) = -( (0.5 - x) * (1.0 - 4.0*z**2) )
      d_basis_3D_dy ( 18, i) = -( (0.5 + x) * (1.0 - 4.0*z**2) )
      d_basis_3D_dy ( 19, i) =  ( (0.5 + x) * (1.0 - 4.0*z**2) )
      d_basis_3D_dy ( 20, i) =  ( (0.5 - x) * (1.0 - 4.0*z**2) )


! --------------------
! z - derivative of Serendipity Basis functions for a 20 Node Cube.

! ----------
! Vertex Nodes - lower z-plane
      d_basis_3D_dz (  1, i) =  ( (-1.0+ 2.0*x) * (-1.0+ 2.0*y)
     1                  * ( 1.0+ 2.0*x + 2.0*y + 4.0*z) ) * 0.25

      d_basis_3D_dz (  2, i) =  ( ( 1.0+ 2.0*x) * (-1.0+ 2.0*y)
     1                  * (-1.0+ 2.0*x - 2.0*y - 4.0*z) ) * 0.25

      d_basis_3D_dz (  3, i) = -( ( 1.0+ 2.0*x) * ( 1.0+ 2.0*y)
     1                  * (-1.0+ 2.0*x + 2.0*y - 4.0*z) ) * 0.25

      d_basis_3D_dz (  4, i) = -( (-1.0+ 2.0*x) * ( 1.0+ 2.0*y)
     1                  * ( 1.0+ 2.0*x - 2.0*y + 4.0*z) ) * 0.25

! ----------
! Vertex Nodes - upper z-plane
      d_basis_3D_dz (  5, i) = -( (-1.0+ 2.0*x) * (-1.0+ 2.0*y)
     1                  * ( 1.0+ 2.0*x + 2.0*y - 4.0*z) ) * 0.25

      d_basis_3D_dz (  6, i) = -( ( 1.0+ 2.0*x) * (-1.0+ 2.0*y)
     1                  * (-1.0+ 2.0*x - 2.0*y + 4.0*z) ) * 0.25

      d_basis_3D_dz (  7, i) =  ( ( 1.0+ 2.0*x) * ( 1.0+ 2.0*y)
     1                  * (-1.0+ 2.0*x + 2.0*y + 4.0*z) ) * 0.25

      d_basis_3D_dz (  8, i) =  ( (-1.0+ 2.0*x) * ( 1.0+ 2.0*y)
     1                  * ( 1.0+ 2.0*x - 2.0*y - 4.0*z) ) * 0.25

! ----------
! x = 0 Nodes
      d_basis_3D_dz (  9, i) = -( (1.0 - 4.0*x**2) * (0.5 - y) )
      d_basis_3D_dz ( 10, i) = -( (1.0 - 4.0*x**2) * (0.5 + y) )
      d_basis_3D_dz ( 11, i) =  ( (1.0 - 4.0*x**2) * (0.5 + y) )
      d_basis_3D_dz ( 12, i) =  ( (1.0 - 4.0*x**2) * (0.5 - y) )

! ----------
! y = 0 Nodes
      d_basis_3D_dz ( 13, i) = -( (0.5 - x) * (1.0 - 4.0*y**2) )
      d_basis_3D_dz ( 14, i) =  ( (0.5 - x) * (1.0 - 4.0*y**2) )
      d_basis_3D_dz ( 15, i) =  ( (0.5 + x) * (1.0 - 4.0*y**2) )
      d_basis_3D_dz ( 16, i) = -( (0.5 + x) * (1.0 - 4.0*y**2) )

! ----------
! z = 0 Nodes
      d_basis_3D_dz ( 17, i) = -8.0 * (0.5 - x)  *  (0.5 - y) * z
      d_basis_3D_dz ( 18, i) = -8.0 * (0.5 + x)  *  (0.5 - y) * z
      d_basis_3D_dz ( 19, i) = -8.0 * (0.5 + x)  *  (0.5 + y) * z
      d_basis_3D_dz ( 20, i) = -8.0 * (0.5 - x)  *  (0.5 + y) * z

! --------------------

          END DO

        END DO

      END DO

! ----------

!$acc data copyin(                                                      &
!$acc&   basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz, wt)

! ----------------------------------------------------------------------
! STEP 2.a) Using "vert_coords_patch" as input, along with
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz" as input, we build
! "v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch" in the volume.
! This involves a large number of small matrix-multiplies, which can be
! done very efficiently on a GPU.

! By the end of this loop nest, if all the quadrature points within the
! volume "v_quadwt_patch" are added up, we should get a good approximation
! for the volume of that zone.


!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    radius1, radius2, tempaa, vert,                               &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, Det)          &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz,        &
!$OMP&    wt, smallnum,                                                 &

!$OMP&    vert_coords_patch,                                            &
!$OMP&    v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch)         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if( do_on_gpu) vector_length(64) present (               &
!$acc&    basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz,        &
!$acc&    wt,                                                           &
!$acc&    vert_coords_patch,                                            &
!$acc&    v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch)         &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    radius1, radius2, tempaa, vert,                               &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, Det)


      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------

           vert ( 1, :) = vert_coords_patch ( ixxm1, iyym1,izzm1,:)
           vert ( 2, :) = vert_coords_patch ( ixx, iyym1, izzm1, :)
           vert ( 3, :) = vert_coords_patch ( ixx, iyy, izzm1, :)
           vert ( 4, :) = vert_coords_patch ( ixxm1, iyy, izzm1, :)

           vert ( 5, :) = vert_coords_patch ( ixxm1, iyym1, izz, :)
           vert ( 6, :) = vert_coords_patch ( ixx, iyym1, izz, :)
           vert ( 7, :) = vert_coords_patch ( ixx, iyy, izz, :)
           vert ( 8, :) = vert_coords_patch ( ixxm1, iyy, izz, :)

          radius1 = SQRT ( vert ( 1, 1)**2
     1                   + vert ( 1, 2)**2  ! Radius of inner face.
     1                   + vert ( 1, 3)**2 )

          radius2 = SQRT ( vert ( 5, 1)**2
     1                   + vert ( 5, 2)**2  ! Radius of outer face.
     1                   + vert ( 5, 3)**2 )

! ------------------------------
! Build the midpoints of the four arcs that make up the x-midplane.

           vert ( 9, :) = 0.5 * ( vert ( 1, :) +  vert ( 2, :))
          tempaa = AMAX1 ( smallnum, SQRT( vert ( 9, 1)**2
     1                                   + vert ( 9, 2)**2
     1                                   + vert ( 9, 3)**2))
          tempaa = radius1 / tempaa
           vert ( 9, :) =  vert ( 9, :) * tempaa

           vert ( 10, :) = 0.5 * ( vert ( 3, :) +  vert ( 4, :))
          tempaa= AMAX1 ( smallnum, SQRT( vert ( 10, 1)**2
     1                                  + vert ( 10, 2)**2
     1                                  + vert ( 10, 3)**2))
          tempaa = radius1 / tempaa
           vert ( 10, :) =  vert ( 10, :) * tempaa

           vert ( 11, :) = 0.5 * ( vert ( 7, :) + vert ( 8, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 11, 1)**2
     1                                  + vert ( 11, 2)**2
     1                                  + vert ( 11, 3)**2))
          tempaa = radius2 / tempaa
           vert ( 11, :) =  vert ( 11, :) * tempaa

           vert ( 12, :) = 0.5 * ( vert ( 5, :) + vert ( 6, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 12, 1)**2
     1                                  + vert ( 12, 2)**2
     1                                  + vert ( 12, 3)**2))
          tempaa = radius2 / tempaa
           vert ( 12, :) =  vert ( 12, :) * tempaa

! ------------------------------
! Build the midpoints of the four arcs that make up the y-midplane.

           vert ( 13, :) = 0.5 * ( vert ( 1, :) +  vert ( 4, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 13, 1)**2
     1                                  + vert ( 13, 2)**2
     1                                  + vert ( 13, 3)**2))
          tempaa = radius1 / tempaa
           vert ( 13, :) =  vert ( 13, :) * tempaa

           vert ( 16, :) = 0.5 * ( vert ( 2, :) +  vert ( 3, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 16, 1)**2
     1                                  + vert ( 16, 2)**2
     1                                  + vert ( 16, 3)**2))
          tempaa = radius1 / tempaa
           vert ( 16, :) =  vert ( 16, :) * tempaa

           vert ( 14, :) = 0.5 * ( vert ( 5, :) +  vert ( 8, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 14, 1)**2
     1                                  + vert ( 14, 2)**2
     1                                  + vert ( 14, 3)**2))
          tempaa = radius2 / tempaa
           vert ( 14, :) =  vert ( 14, :) * tempaa

           vert ( 15, :) = 0.5 * ( vert ( 6, :) +  vert ( 7, :))
          tempaa =AMAX1 ( smallnum, SQRT( vert ( 15, 1)**2
     1                                  + vert ( 15, 2)**2
     1                                  + vert ( 15, 3)**2))
          tempaa = radius2 / tempaa
           vert ( 15, :) =  vert ( 15, :) * tempaa

! ------------------------------
! Build the midpoints of the four straight edges that make up the z-midplane.

           vert ( 17, :) = 0.5 * ( vert ( 1, :) +  vert ( 5, :))

           vert ( 18, :) = 0.5 * ( vert ( 2, :) +  vert ( 6, :))

           vert ( 19, :) = 0.5 * ( vert ( 3, :) +  vert ( 7, :))

           vert ( 20, :) = 0.5 * ( vert ( 4, :) +  vert ( 8, :))

! ------------------------------
! Now we extract x,y,z-components from "vert" and do the MATMUL with
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz".

! We are already in a position to load in the three components of the
! nodal positions in this volume "v_quadpts_patch ( ixx, iyy, izz, 1:3, :)".

! We also build the vector ( Ax, Ay, Az) which is the partial derivative
! of the position vector w.r.t. the x-coordinate in the reference element.

! We also build the vector ( Bx, By, Bz) which is the partial derivative
! of the position vector w.r.t. the y-coordinate in the reference element.

! We also build the vector ( Cx, Cy, Cz) which is the partial derivative
! of the position vector w.r.t. the z-coordinate in the reference element.

! ----------
! Extract all the global x-coordinates and their x, y and z derivatives w.r.t.
! the reference element by multiplying "vert ( 1:20, 1)" with
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz".
! Store in "v_quadpts_patch ( ixx, iyy, izz, 1, :)".

          v_quadpts_patch ( ixx, iyy, izz, 1, :) =
     1       vert (  1, 1) * basis_3D (  1, :)
     1     + vert (  2, 1) * basis_3D (  2, :)
     1     + vert (  3, 1) * basis_3D (  3, :)
     1     + vert (  4, 1) * basis_3D (  4, :)
     1     + vert (  5, 1) * basis_3D (  5, :)
     1     + vert (  6, 1) * basis_3D (  6, :)
     1     + vert (  7, 1) * basis_3D (  7, :)
     1     + vert (  8, 1) * basis_3D (  8, :)
     1     + vert (  9, 1) * basis_3D (  9, :)
     1     + vert ( 10, 1) * basis_3D ( 10, :)
     1     + vert ( 11, 1) * basis_3D ( 11, :)
     1     + vert ( 12, 1) * basis_3D ( 12, :)
     1     + vert ( 13, 1) * basis_3D ( 13, :)
     1     + vert ( 14, 1) * basis_3D ( 14, :)
     1     + vert ( 15, 1) * basis_3D ( 15, :)
     1     + vert ( 16, 1) * basis_3D ( 16, :)
     1     + vert ( 17, 1) * basis_3D ( 17, :)
     1     + vert ( 18, 1) * basis_3D ( 18, :)
     1     + vert ( 19, 1) * basis_3D ( 19, :)
     1     + vert ( 20, 1) * basis_3D ( 20, :)

          Ax ( :) =
     1       vert (  1, 1) * d_basis_3D_dx (  1, :)
     1     + vert (  2, 1) * d_basis_3D_dx (  2, :)
     1     + vert (  3, 1) * d_basis_3D_dx (  3, :)
     1     + vert (  4, 1) * d_basis_3D_dx (  4, :)
     1     + vert (  5, 1) * d_basis_3D_dx (  5, :)
     1     + vert (  6, 1) * d_basis_3D_dx (  6, :)
     1     + vert (  7, 1) * d_basis_3D_dx (  7, :)
     1     + vert (  8, 1) * d_basis_3D_dx (  8, :)
     1     + vert (  9, 1) * d_basis_3D_dx (  9, :)
     1     + vert ( 10, 1) * d_basis_3D_dx ( 10, :)
     1     + vert ( 11, 1) * d_basis_3D_dx ( 11, :)
     1     + vert ( 12, 1) * d_basis_3D_dx ( 12, :)
     1     + vert ( 13, 1) * d_basis_3D_dx ( 13, :)
     1     + vert ( 14, 1) * d_basis_3D_dx ( 14, :)
     1     + vert ( 15, 1) * d_basis_3D_dx ( 15, :)
     1     + vert ( 16, 1) * d_basis_3D_dx ( 16, :)
     1     + vert ( 17, 1) * d_basis_3D_dx ( 17, :)
     1     + vert ( 18, 1) * d_basis_3D_dx ( 18, :)
     1     + vert ( 19, 1) * d_basis_3D_dx ( 19, :)
     1     + vert ( 20, 1) * d_basis_3D_dx ( 20, :)

          Bx ( :) =
     1       vert (  1, 1) * d_basis_3D_dy (  1, :)
     1     + vert (  2, 1) * d_basis_3D_dy (  2, :)
     1     + vert (  3, 1) * d_basis_3D_dy (  3, :)
     1     + vert (  4, 1) * d_basis_3D_dy (  4, :)
     1     + vert (  5, 1) * d_basis_3D_dy (  5, :)
     1     + vert (  6, 1) * d_basis_3D_dy (  6, :)
     1     + vert (  7, 1) * d_basis_3D_dy (  7, :)
     1     + vert (  8, 1) * d_basis_3D_dy (  8, :)
     1     + vert (  9, 1) * d_basis_3D_dy (  9, :)
     1     + vert ( 10, 1) * d_basis_3D_dy ( 10, :)
     1     + vert ( 11, 1) * d_basis_3D_dy ( 11, :)
     1     + vert ( 12, 1) * d_basis_3D_dy ( 12, :)
     1     + vert ( 13, 1) * d_basis_3D_dy ( 13, :)
     1     + vert ( 14, 1) * d_basis_3D_dy ( 14, :)
     1     + vert ( 15, 1) * d_basis_3D_dy ( 15, :)
     1     + vert ( 16, 1) * d_basis_3D_dy ( 16, :)
     1     + vert ( 17, 1) * d_basis_3D_dy ( 17, :)
     1     + vert ( 18, 1) * d_basis_3D_dy ( 18, :)
     1     + vert ( 19, 1) * d_basis_3D_dy ( 19, :)
     1     + vert ( 20, 1) * d_basis_3D_dy ( 20, :)

          Cx ( :) =
     1       vert (  1, 1) * d_basis_3D_dz (  1, :)
     1     + vert (  2, 1) * d_basis_3D_dz (  2, :)
     1     + vert (  3, 1) * d_basis_3D_dz (  3, :)
     1     + vert (  4, 1) * d_basis_3D_dz (  4, :)
     1     + vert (  5, 1) * d_basis_3D_dz (  5, :)
     1     + vert (  6, 1) * d_basis_3D_dz (  6, :)
     1     + vert (  7, 1) * d_basis_3D_dz (  7, :)
     1     + vert (  8, 1) * d_basis_3D_dz (  8, :)
     1     + vert (  9, 1) * d_basis_3D_dz (  9, :)
     1     + vert ( 10, 1) * d_basis_3D_dz ( 10, :)
     1     + vert ( 11, 1) * d_basis_3D_dz ( 11, :)
     1     + vert ( 12, 1) * d_basis_3D_dz ( 12, :)
     1     + vert ( 13, 1) * d_basis_3D_dz ( 13, :)
     1     + vert ( 14, 1) * d_basis_3D_dz ( 14, :)
     1     + vert ( 15, 1) * d_basis_3D_dz ( 15, :)
     1     + vert ( 16, 1) * d_basis_3D_dz ( 16, :)
     1     + vert ( 17, 1) * d_basis_3D_dz ( 17, :)
     1     + vert ( 18, 1) * d_basis_3D_dz ( 18, :)
     1     + vert ( 19, 1) * d_basis_3D_dz ( 19, :)
     1     + vert ( 20, 1) * d_basis_3D_dz ( 20, :)

! ----------
! Extract all the global y-coordinates and their x, y and z derivatives w.r.t.
! the reference element by multiplying "vert ( 1:20, 2)" with
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz".
! Store in "v_quadpts_patch ( ixx, iyy, izz, 2, :)".

          v_quadpts_patch ( ixx, iyy, izz, 2, :) =
     1       vert (  1, 2) * basis_3D (  1, :)
     1     + vert (  2, 2) * basis_3D (  2, :)
     1     + vert (  3, 2) * basis_3D (  3, :)
     1     + vert (  4, 2) * basis_3D (  4, :)
     1     + vert (  5, 2) * basis_3D (  5, :)
     1     + vert (  6, 2) * basis_3D (  6, :)
     1     + vert (  7, 2) * basis_3D (  7, :)
     1     + vert (  8, 2) * basis_3D (  8, :)
     1     + vert (  9, 2) * basis_3D (  9, :)
     1     + vert ( 10, 2) * basis_3D ( 10, :)
     1     + vert ( 11, 2) * basis_3D ( 11, :)
     1     + vert ( 12, 2) * basis_3D ( 12, :)
     1     + vert ( 13, 2) * basis_3D ( 13, :)
     1     + vert ( 14, 2) * basis_3D ( 14, :)
     1     + vert ( 15, 2) * basis_3D ( 15, :)
     1     + vert ( 16, 2) * basis_3D ( 16, :)
     1     + vert ( 17, 2) * basis_3D ( 17, :)
     1     + vert ( 18, 2) * basis_3D ( 18, :)
     1     + vert ( 19, 2) * basis_3D ( 19, :)
     1     + vert ( 20, 2) * basis_3D ( 20, :)

          Ay ( :) =
     1       vert (  1, 2) * d_basis_3D_dx (  1, :)
     1     + vert (  2, 2) * d_basis_3D_dx (  2, :)
     1     + vert (  3, 2) * d_basis_3D_dx (  3, :)
     1     + vert (  4, 2) * d_basis_3D_dx (  4, :)
     1     + vert (  5, 2) * d_basis_3D_dx (  5, :)
     1     + vert (  6, 2) * d_basis_3D_dx (  6, :)
     1     + vert (  7, 2) * d_basis_3D_dx (  7, :)
     1     + vert (  8, 2) * d_basis_3D_dx (  8, :)
     1     + vert (  9, 2) * d_basis_3D_dx (  9, :)
     1     + vert ( 10, 2) * d_basis_3D_dx ( 10, :)
     1     + vert ( 11, 2) * d_basis_3D_dx ( 11, :)
     1     + vert ( 12, 2) * d_basis_3D_dx ( 12, :)
     1     + vert ( 13, 2) * d_basis_3D_dx ( 13, :)
     1     + vert ( 14, 2) * d_basis_3D_dx ( 14, :)
     1     + vert ( 15, 2) * d_basis_3D_dx ( 15, :)
     1     + vert ( 16, 2) * d_basis_3D_dx ( 16, :)
     1     + vert ( 17, 2) * d_basis_3D_dx ( 17, :)
     1     + vert ( 18, 2) * d_basis_3D_dx ( 18, :)
     1     + vert ( 19, 2) * d_basis_3D_dx ( 19, :)
     1     + vert ( 20, 2) * d_basis_3D_dx ( 20, :)

          By ( :) =
     1       vert (  1, 2) * d_basis_3D_dy (  1, :)
     1     + vert (  2, 2) * d_basis_3D_dy (  2, :)
     1     + vert (  3, 2) * d_basis_3D_dy (  3, :)
     1     + vert (  4, 2) * d_basis_3D_dy (  4, :)
     1     + vert (  5, 2) * d_basis_3D_dy (  5, :)
     1     + vert (  6, 2) * d_basis_3D_dy (  6, :)
     1     + vert (  7, 2) * d_basis_3D_dy (  7, :)
     1     + vert (  8, 2) * d_basis_3D_dy (  8, :)
     1     + vert (  9, 2) * d_basis_3D_dy (  9, :)
     1     + vert ( 10, 2) * d_basis_3D_dy ( 10, :)
     1     + vert ( 11, 2) * d_basis_3D_dy ( 11, :)
     1     + vert ( 12, 2) * d_basis_3D_dy ( 12, :)
     1     + vert ( 13, 2) * d_basis_3D_dy ( 13, :)
     1     + vert ( 14, 2) * d_basis_3D_dy ( 14, :)
     1     + vert ( 15, 2) * d_basis_3D_dy ( 15, :)
     1     + vert ( 16, 2) * d_basis_3D_dy ( 16, :)
     1     + vert ( 17, 2) * d_basis_3D_dy ( 17, :)
     1     + vert ( 18, 2) * d_basis_3D_dy ( 18, :)
     1     + vert ( 19, 2) * d_basis_3D_dy ( 19, :)
     1     + vert ( 20, 2) * d_basis_3D_dy ( 20, :)

          Cy ( :) =
     1       vert (  1, 2) * d_basis_3D_dz (  1, :)
     1     + vert (  2, 2) * d_basis_3D_dz (  2, :)
     1     + vert (  3, 2) * d_basis_3D_dz (  3, :)
     1     + vert (  4, 2) * d_basis_3D_dz (  4, :)
     1     + vert (  5, 2) * d_basis_3D_dz (  5, :)
     1     + vert (  6, 2) * d_basis_3D_dz (  6, :)
     1     + vert (  7, 2) * d_basis_3D_dz (  7, :)
     1     + vert (  8, 2) * d_basis_3D_dz (  8, :)
     1     + vert (  9, 2) * d_basis_3D_dz (  9, :)
     1     + vert ( 10, 2) * d_basis_3D_dz ( 10, :)
     1     + vert ( 11, 2) * d_basis_3D_dz ( 11, :)
     1     + vert ( 12, 2) * d_basis_3D_dz ( 12, :)
     1     + vert ( 13, 2) * d_basis_3D_dz ( 13, :)
     1     + vert ( 14, 2) * d_basis_3D_dz ( 14, :)
     1     + vert ( 15, 2) * d_basis_3D_dz ( 15, :)
     1     + vert ( 16, 2) * d_basis_3D_dz ( 16, :)
     1     + vert ( 17, 2) * d_basis_3D_dz ( 17, :)
     1     + vert ( 18, 2) * d_basis_3D_dz ( 18, :)
     1     + vert ( 19, 2) * d_basis_3D_dz ( 19, :)
     1     + vert ( 20, 2) * d_basis_3D_dz ( 20, :)

! ----------
! Extract all the global z-coordinates and their x, y and z derivatives w.r.t.
! the reference element by multiplying "vert ( 1:20, 3)" with
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz".
! Store in "v_quadpts_patch ( ixx, iyy, izz, 3, :)".

          v_quadpts_patch ( ixx, iyy, izz, 3, :) =
     1       vert (  1, 3) * basis_3D (  1, :)
     1     + vert (  2, 3) * basis_3D (  2, :)
     1     + vert (  3, 3) * basis_3D (  3, :)
     1     + vert (  4, 3) * basis_3D (  4, :)
     1     + vert (  5, 3) * basis_3D (  5, :)
     1     + vert (  6, 3) * basis_3D (  6, :)
     1     + vert (  7, 3) * basis_3D (  7, :)
     1     + vert (  8, 3) * basis_3D (  8, :)
     1     + vert (  9, 3) * basis_3D (  9, :)
     1     + vert ( 10, 3) * basis_3D ( 10, :)
     1     + vert ( 11, 3) * basis_3D ( 11, :)
     1     + vert ( 12, 3) * basis_3D ( 12, :)
     1     + vert ( 13, 3) * basis_3D ( 13, :)
     1     + vert ( 14, 3) * basis_3D ( 14, :)
     1     + vert ( 15, 3) * basis_3D ( 15, :)
     1     + vert ( 16, 3) * basis_3D ( 16, :)
     1     + vert ( 17, 3) * basis_3D ( 17, :)
     1     + vert ( 18, 3) * basis_3D ( 18, :)
     1     + vert ( 19, 3) * basis_3D ( 19, :)
     1     + vert ( 20, 3) * basis_3D ( 20, :)

          Az ( :) =
     1       vert (  1, 3) * d_basis_3D_dx (  1, :)
     1     + vert (  2, 3) * d_basis_3D_dx (  2, :)
     1     + vert (  3, 3) * d_basis_3D_dx (  3, :)
     1     + vert (  4, 3) * d_basis_3D_dx (  4, :)
     1     + vert (  5, 3) * d_basis_3D_dx (  5, :)
     1     + vert (  6, 3) * d_basis_3D_dx (  6, :)
     1     + vert (  7, 3) * d_basis_3D_dx (  7, :)
     1     + vert (  8, 3) * d_basis_3D_dx (  8, :)
     1     + vert (  9, 3) * d_basis_3D_dx (  9, :)
     1     + vert ( 10, 3) * d_basis_3D_dx ( 10, :)
     1     + vert ( 11, 3) * d_basis_3D_dx ( 11, :)
     1     + vert ( 12, 3) * d_basis_3D_dx ( 12, :)
     1     + vert ( 13, 3) * d_basis_3D_dx ( 13, :)
     1     + vert ( 14, 3) * d_basis_3D_dx ( 14, :)
     1     + vert ( 15, 3) * d_basis_3D_dx ( 15, :)
     1     + vert ( 16, 3) * d_basis_3D_dx ( 16, :)
     1     + vert ( 17, 3) * d_basis_3D_dx ( 17, :)
     1     + vert ( 18, 3) * d_basis_3D_dx ( 18, :)
     1     + vert ( 19, 3) * d_basis_3D_dx ( 19, :)
     1     + vert ( 20, 3) * d_basis_3D_dx ( 20, :)

          Bz ( :) =
     1       vert (  1, 3) * d_basis_3D_dy (  1, :)
     1     + vert (  2, 3) * d_basis_3D_dy (  2, :)
     1     + vert (  3, 3) * d_basis_3D_dy (  3, :)
     1     + vert (  4, 3) * d_basis_3D_dy (  4, :)
     1     + vert (  5, 3) * d_basis_3D_dy (  5, :)
     1     + vert (  6, 3) * d_basis_3D_dy (  6, :)
     1     + vert (  7, 3) * d_basis_3D_dy (  7, :)
     1     + vert (  8, 3) * d_basis_3D_dy (  8, :)
     1     + vert (  9, 3) * d_basis_3D_dy (  9, :)
     1     + vert ( 10, 3) * d_basis_3D_dy ( 10, :)
     1     + vert ( 11, 3) * d_basis_3D_dy ( 11, :)
     1     + vert ( 12, 3) * d_basis_3D_dy ( 12, :)
     1     + vert ( 13, 3) * d_basis_3D_dy ( 13, :)
     1     + vert ( 14, 3) * d_basis_3D_dy ( 14, :)
     1     + vert ( 15, 3) * d_basis_3D_dy ( 15, :)
     1     + vert ( 16, 3) * d_basis_3D_dy ( 16, :)
     1     + vert ( 17, 3) * d_basis_3D_dy ( 17, :)
     1     + vert ( 18, 3) * d_basis_3D_dy ( 18, :)
     1     + vert ( 19, 3) * d_basis_3D_dy ( 19, :)
     1     + vert ( 20, 3) * d_basis_3D_dy ( 20, :)

          Cz ( :) =
     1       vert (  1, 3) * d_basis_3D_dz (  1, :)
     1     + vert (  2, 3) * d_basis_3D_dz (  2, :)
     1     + vert (  3, 3) * d_basis_3D_dz (  3, :)
     1     + vert (  4, 3) * d_basis_3D_dz (  4, :)
     1     + vert (  5, 3) * d_basis_3D_dz (  5, :)
     1     + vert (  6, 3) * d_basis_3D_dz (  6, :)
     1     + vert (  7, 3) * d_basis_3D_dz (  7, :)
     1     + vert (  8, 3) * d_basis_3D_dz (  8, :)
     1     + vert (  9, 3) * d_basis_3D_dz (  9, :)
     1     + vert ( 10, 3) * d_basis_3D_dz ( 10, :)
     1     + vert ( 11, 3) * d_basis_3D_dz ( 11, :)
     1     + vert ( 12, 3) * d_basis_3D_dz ( 12, :)
     1     + vert ( 13, 3) * d_basis_3D_dz ( 13, :)
     1     + vert ( 14, 3) * d_basis_3D_dz ( 14, :)
     1     + vert ( 15, 3) * d_basis_3D_dz ( 15, :)
     1     + vert ( 16, 3) * d_basis_3D_dz ( 16, :)
     1     + vert ( 17, 3) * d_basis_3D_dz ( 17, :)
     1     + vert ( 18, 3) * d_basis_3D_dz ( 18, :)
     1     + vert ( 19, 3) * d_basis_3D_dz ( 19, :)
     1     + vert ( 20, 3) * d_basis_3D_dz ( 20, :)

! ------------------------------
! Now we obtain the Determinant of "C dot ( A cross B)" and store in "Det". We
! take its absolute value. Along with the weight "wt",
! that will give us the nodal weights "v_quadwt_patch". We also save "A, B, C"
! vectors in "v_quadtangent_patch ( ixx, iyy, izz, 1:3, 1:3, :)"; useful for ADER.

          Dx ( :) = Ay ( :) * Bz ( :) - Az ( :) * By ( :)
          Dy ( :) = Az ( :) * Bx ( :) - Ax ( :) * Bz ( :)
          Dz ( :) = Ax ( :) * By ( :) - Ay ( :) * Bx ( :)

          Det ( :) = ABS ( Cx ( :) * Dx ( :)
     1                   + Cy ( :) * Dy ( :)
     1                   + Cz ( :) * Dz ( :) )

          v_quadwt_patch ( ixx, iyy, izz, :) = wt ( :) * Det ( :)

! Derivative of position vector w.r.t. x-coordinate in the reference element.
          v_quadtangent_patch ( ixx, iyy, izz, 1, 1, :) = Ax ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 2, 1, :) = Ay ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 3, 1, :) = Az ( :)

! Derivative of position vector w.r.t. y-coordinate in the reference element.
          v_quadtangent_patch ( ixx, iyy, izz, 1, 2, :) = Bx ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 2, 2, :) = By ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 3, 2, :) = Bz ( :)

! Derivative of position vector w.r.t. z-coordinate in the reference element.
          v_quadtangent_patch ( ixx, iyy, izz, 1, 3, :) = Cx ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 2, 3, :) = Cy ( :)
          v_quadtangent_patch ( ixx, iyy, izz, 3, 3, :) = Cz ( :)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! STEP 2.b) Normalize "v_quadwt_patch" weights to add up to unity.
! Also obtain "v_volume_patch".

!$acc kernels if( do_on_gpu)
      v_volume_patch ( :, :, :) = 0.0
!$acc end kernels

      DO i = 1, NUMORDER * NUMORDER * NUMORDER

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, i,              &
!$OMP&    v_volume_patch, v_quadwt_patch)                               &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if( do_on_gpu) vector_length(64) present (               &
!$acc&    v_volume_patch, v_quadwt_patch)                               &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

! ------------------------------

          v_volume_patch ( ixx, iyy, izz) =
     1    v_volume_patch ( ixx, iyy, izz)
     1  + v_quadwt_patch ( ixx, iyy, izz, i)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO i = 1, NUMORDER * NUMORDER" loop.

! --------------------------------------------------
! STEP 2.c) Make sure "v_quadwt_patch" weights add up to unity.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    v_volume_patch, v_quadwt_patch)                               &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if( do_on_gpu) vector_length(64) present (               &
!$acc&    v_volume_patch, v_quadwt_patch)                               &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

! ------------------------------

          v_volume_patch ( ixx, iyy, izz) = AMAX1 ( smallnum,
     1                   v_volume_patch ( ixx, iyy, izz) )

          v_quadwt_patch ( ixx, iyy, izz, :) =
     1    v_quadwt_patch ( ixx, iyy, izz, :)
     1  / v_volume_patch ( ixx, iyy, izz)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! STEP 3) Obtain the centroid "v_centroid_patch" for each zone.
! "v_centroid_temp" is used as a large temporary variable.

!$acc kernels if ( do_on_gpu)
      v_centroid_patch ( :, :, :, :) = 0.0
!$acc end kernels

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    vert)                                                         &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    vert_coords_patch, v_centroid_patch,                          &
!$OMP&    v_quadwt_patch, v_quadpts_patch)                              &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if( do_on_gpu) vector_length(64) present (               &
!$acc&    vert_coords_patch, v_centroid_patch,                          &
!$acc&    v_quadwt_patch, v_quadpts_patch)                              &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    vert)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------

           vert ( 1, :) = vert_coords_patch ( ixxm1, iyym1,izzm1,:)
           vert ( 2, :) = vert_coords_patch ( ixx, iyym1, izzm1, :)
           vert ( 3, :) = vert_coords_patch ( ixx, iyy, izzm1, :)
           vert ( 4, :) = vert_coords_patch ( ixxm1, iyy, izzm1, :)

           vert ( 5, :) = vert_coords_patch ( ixxm1, iyym1, izz, :)
           vert ( 6, :) = vert_coords_patch ( ixx, iyym1, izz, :)
           vert ( 7, :) = vert_coords_patch ( ixx, iyy, izz, :)
           vert ( 8, :) = vert_coords_patch ( ixxm1, iyy, izz, :)

! ------------------------------
! Build rough approximation of the centroid by averaging the 8 vertices.

          v_centroid_patch ( ixx, iyy, izz, 1) = 0.125 * (
     1     vert ( 1, 1) +  vert ( 2, 1) +  vert ( 3, 1)
     1  +  vert ( 4, 1) +  vert ( 5, 1) +  vert ( 6, 1)
     1  +  vert ( 7, 1) +  vert ( 8, 1) )

          v_centroid_patch ( ixx, iyy, izz, 2) = 0.125 * (
     1     vert ( 1, 2) +  vert ( 2, 2) +  vert ( 3, 2)
     1  +  vert ( 4, 2) +  vert ( 5, 2) +  vert ( 6, 2)
     1  +  vert ( 7, 2) +  vert ( 8, 2) )

          v_centroid_patch ( ixx, iyy, izz, 3) = 0.125 * (
     1     vert ( 1, 3) +  vert ( 2, 3) +  vert ( 3, 3)
     1  +  vert ( 4, 3) +  vert ( 5, 3) +  vert ( 6, 3)
     1  +  vert ( 7, 3) +  vert ( 8, 3) )

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO
! ------------------------------
! Obtain "v_centroid_temp".

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&                                    v_centroid_temp)              &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)
      DO izz = 1 - NGhost, nr_zones + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
          DO ixx = 1 - NGhost, ntp_zones + NGhost
          v_centroid_temp ( ixx, iyy, izz, :) = 0.0
          END DO
        END DO
      END DO
!$acc end parallel

! ----------

      DO i = 1, NUMORDER * NUMORDER * NUMORDER

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, i, v_centroid_temp,     &
!$OMP&             v_quadwt_patch, v_quadpts_patch, v_centroid_patch)   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&             v_centroid_temp,                                     &
!$acc&             v_quadwt_patch, v_quadpts_patch, v_centroid_patch)   &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! ------------------------------
! x-coord wrt "centroid".
          v_centroid_temp ( ixx, iyy, izz, 1) =
     1                     v_centroid_temp ( ixx, iyy, izz, 1)
     1                    + v_quadwt_patch ( ixx, iyy, izz, i)
     1              * ( v_quadpts_patch ( ixx, iyy, izz, 1, i)
     1                - v_centroid_patch ( ixx, iyy, izz, 1))

! y-coord wrt "centroid".
          v_centroid_temp ( ixx, iyy, izz, 2) =
     1                     v_centroid_temp ( ixx, iyy, izz, 2)
     1                    + v_quadwt_patch ( ixx, iyy, izz, i)
     1              * ( v_quadpts_patch ( ixx, iyy, izz, 2, i)
     1                - v_centroid_patch ( ixx, iyy, izz, 2))

! z-coord wrt "centroid".
          v_centroid_temp ( ixx, iyy, izz, 3) =
     1                     v_centroid_temp ( ixx, iyy, izz, 3)
     1                    + v_quadwt_patch ( ixx, iyy, izz, i)
     1              * ( v_quadpts_patch ( ixx, iyy, izz, 3, i)
     1                - v_centroid_patch ( ixx, iyy, izz, 3))

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO i = 1, NUMORDER * NUMORDER" loop.

! ------------------------------
! Shift the conjectured "v_centroid_patch" by "v_centroid_temp" to get the precise
! centroid.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, v_centroid_temp,        &
!$OMP&             v_centroid_patch)                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&                     v_centroid_temp, v_centroid_patch)           &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------
          v_centroid_patch ( ixx, iyy, izz, :) =
     1    v_centroid_patch ( ixx, iyy, izz, :)
     1  + v_centroid_temp  ( ixx, iyy, izz, :)

! ----------
          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.
        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.
!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------

#if (DEBUG_GEOM == 1)

! ------------------------------

!$acc update host (                                                     &
!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch, v_centroid_temp)

! --------------------------------------------------

! We should check to make sure that the centroid lies within the zone.
! This inexpensive test is done serially, for now.

! Commented out because on a GPU, this test would be a deadly bottleneck.

      DO ixx = 1 - NGhost, ntp_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixxm1 = ixx - 1
        iyym1 = iyy - 1

          i = 0

          DO izz = 1 - NGhost, nr_zones + NGhost

          izzm1 = izz - 1

! ------------------------------

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

           vert ( 1, :) = vert_coords_patch ( 1, 1, izzm1, :)

           vert ( 5, :) = vert_coords_patch ( 1, 1, izz, :)

          radius1 = SQRT (  vert ( 1, 1)**2
     1                   +  vert ( 1, 2)**2  ! Radius of inner face.
     1                   +  vert ( 1, 3)**2 )

          radius2 = SQRT (  vert ( 5, 1)**2
     1                   +  vert ( 5, 2)**2  ! Radius of outer face.
     1                   +  vert ( 5, 3)**2 )

          tempaa = SQRT ( v_centroid_patch ( ixx, iyy, izz, 1)**2
     1                  + v_centroid_patch ( ixx, iyy, izz, 2)**2
     1                  + v_centroid_patch ( ixx, iyy, izz, 3)**2 )

          IF ( ( tempaa .LT. radius1)
     1    .OR. ( tempaa .GT. radius2)) THEN
          i = izz
          EXIT
          END IF

          END IF

! ------------------------------

          END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

        IF ( i .NE. 0) THEN
        WRITE ( 6, *)"Stopping in COORDS_TO_WT_QUADPTS_VOLUME_ISO2"
        WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, i
        STOP
        END IF

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

! ------------------------------

#endif

! --------------------------------------------------
! Make several debugging diagnostics, if that is called for.

#if (DEBUG_GEOM == 1)

! ------------------------------
! First test consists of demanding that the volume of any one complete
! sector is one-sixth the volume of the spherical shell.

           vert ( 1, :) = vert_coords_patch ( 1, 1, 0, :)

           vert ( 5, :) = vert_coords_patch ( 1, 1, nr_zones, :)

          radius1 = SQRT (  vert ( 1, 1)**2
     1                   +  vert ( 1, 2)**2  ! Radius of inner face.
     1                   +  vert ( 1, 3)**2 )

          radius2 = SQRT (  vert ( 5, 1)**2
     1                   +  vert ( 5, 2)**2  ! Radius of outer face.
     1                   +  vert ( 5, 3)**2 )

          tempaa = (2.0 * pi / 9.0) * ( radius2**3
     1                                       - radius1**3)  ! Analytic volume.

          tempab = 0.0

          DO ixx = 1, ntp_zones
            DO iyy = 1, ntp_zones
              DO izz = 1, nr_zones
              tempab = tempab
     1                      + v_volume_patch ( ixx, iyy, izz)
              END DO
            END DO
          END DO

          WRITE ( 6, *)" "
          WRITE ( 6, *)"Analytic volume, numeric volume, Diff = ",
     1                 tempaa, tempab,
     1           ABS ( tempaa - tempab)
          WRITE ( 6, *)" "

! ------------------------------
! Since the evaluation of the centroid is a linear process, evaluating the
! centroid all over again should yield zero shift in the centroid's location.
! This is our second test.


      DO ixx = 1 - NGhost, ntp_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixxm1 = ixx - 1
        iyym1 = iyy - 1

          DO izz = 1 - NGhost, nr_zones + NGhost

          izzm1 = izz - 1

! ------------------------------

          qwt ( 1, :) = v_quadwt_patch ( ixx, iyy, izz, :)  ! Vector of weights.

          vec_out ( :) = v_quadpts_patch ( ixx, iyy, izz, 1, :)  ! x-coord wrt "centroid".
     1                 - v_centroid_patch ( ixx, iyy, izz, 1)

          q11 = MATMUL ( qwt, vec_out)
          v_centroid_temp ( ixx, iyy, izz, 1) = q11 ( 1)

          vec_out ( :) = v_quadpts_patch ( ixx, iyy, izz, 2, :)  ! y-coord wrt "centroid".
     1                 - v_centroid_patch ( ixx, iyy, izz, 2)

          q11 = MATMUL ( qwt, vec_out)
          v_centroid_temp ( ixx, iyy, izz, 2) = q11 ( 1)

          vec_out ( :) = v_quadpts_patch ( ixx, iyy, izz, 3, :)  ! z-coord wrt "centroid".
     1                 - v_centroid_patch ( ixx, iyy, izz, 3)

          q11 = MATMUL ( qwt, vec_out)
          v_centroid_temp ( ixx, iyy, izz, 3) = q11 ( 1)

! ------------------------------

          END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.


! ------------------------------

      tempaa = 0.0

      DO ixx = 1 - NGhost, ntp_zones + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
          DO izz = 1 - NGhost, nr_zones + NGhost

          tempaa = tempaa
     1           + SQRT ( v_centroid_temp ( ixx, iyy, izz, 1)**2
     1                  + v_centroid_temp ( ixx, iyy, izz, 2)**2
     1                  + v_centroid_temp ( ixx, iyy, izz, 3)**2 )

          END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.
        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

      WRITE ( 6, *)" "
      WRITE ( 6, *)"On second round, shift in centroid = ",
     1      tempaa / REAL ( nr_zones * ntp_zones * ntp_zones)
      WRITE ( 6, *)"Must be very close to machine zero."
      WRITE ( 6, *)" "

! ------------------------------

#endif

! Each "data create" needs a matching "end data".

!$acc end data

!$acc end data

! This is because we are releasing the allocation on :
! "basis_3D, d_basis_3D_dx, d_basis_3D_dy, d_basis_3D_dz, wt".

! ----------------------------------------------------------------------

      END SUBROUTINE COORDS_TO_WT_QUADPTS_VOLUME_ISO2


