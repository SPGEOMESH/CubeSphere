






      SUBROUTINE COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara and Sethu Subramanyan

! This subroutine takes in "vert_coords_patch" which contains the coordinates
! of the vertices of one particular patch of the cube sphere mesh. It also takes in
! "imask_active_patch" which is 1 for good zones, 0 for bad ones. It then performs
! isoparametric mapping of degree two in each of the three faces of the mesh.

! This version is specialized for reconstructing the normal components in the
! faces of the mesh. Therefore, even the dimensioning of the "xb_*, yb_*, zb_*"
! variables is different from the "xf_*, yf_*, zf_*" variables that are initialized
! in "COORDS_TO_WT_QUADPTS_FACES_ISO2".

! It returns the locations of the quadrature points in the faces:-
! "xb_quadpts_patch, yb_quadpts_patch, zb_quadpts_patch".

! It also returns the weights associated with the quadrature points in the faces:-
! "xb_quadwt_patch, yb_quadwt_patch, zb_quadwt_patch".

! It also returns the unit normal vectors to the faces. Only 1 such vector is
! needed for the x- and y-faces, because they are flat. But the z-face
! is curved and needs the full set of unit normal vectors.

! The location of quadrature points in each of the faces is as follows:-

!---------------------------------------------------------------------
!                 Schematic of Quadrature Points inside a Cube
!
!                          4 •-----------------------------• 3
!                           /|           q4               /|
!                          / |                           / |
!                         /  |                          /  |
!                        /q3 |         q1            q2/   |
!                       /    |                        /q4  |
!                      / 1   |       q5            2 /     |
!                   2 •----------------------------•4    q2|
!                     |      |        q2         3 |       |
!                     |      |                     |   q1  |
!                     |      |                     |       |
!                     |      •---------------------|-------• 2
!                     |q5   /        q1          q4|q3    /
!                     |    /                       |   q5/
!                     |   /                        |    /
!                     |  /                         |   /
!                     | /           q3             |  /
!                     |/___________________________| /
!                    1 •------------------------- 4 • 1
!
!
!
!
!
!                                         z   y
!                                         ↑  /
!                                         | /
!                                         |/
!                                         •-------------> x
!
!
!
! Here q1...5 are the facial quadrature points. Where, q1 is the central one
! with weight 1.0/3.0 and q2, q3, q4, q5 are the other quadrature points with
! weight 1.0/6.0.

! In this schematic only four vertices (1, 2, 3, 4) of each faces are shown for simplicity.
! But Vertex 5 is between vertices 1 and 2.
! Vertex 6 is between vertices 2 and 3. Vertex 7 is between 3 and 4.
! Vertex 8 is between 4 and 1.
!------------------------------------------------------------
!
!                                          7
!                           4•-----------------------------•3
!                           /            q4               /
!                          /                             /
!                        8/                             /    ---> z-face
!                        /q3           q1            q2/ 6
!                       /                             /
!                      /             q5              /
!                    1•----------------------------•2
!                                   5
!------------------------------------------------------------
!                                    6
!                   2 •----------------------------•3
!                     |              q2            |
!                     |                            |
!                     |                            |
!                     |q5                        q4|
!                    5|              q1            |7    ---> y-face
!                     |                            |
!                     |                            |
!                     |                            |
!                     |              q3            |
!                     |____________________________|
!                    1              8              4
!
!------------------------------------------------------------
!
!                                                • 3
!                                               /|
!                                              / |
!                                             /  |
!                                           7/   |
!                                           /q4  |
!                                          /     |6
!                                        •4    q2|       ---> x-face
!                                        |       |
!                                        |   q1  |
!                                        |       |
!                                       8|       • 2
!                                        |q3    /
!                                        |   q5/
!                                        |    /5
!                                        |   /
!                                        |  /
!                                        | /
!                                         • 1
!
!------------------------------------------------------------
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
! "NBFLD_QUAD_PTS" number of quadrature points in the x-face of a zone.
      REAL xb_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the x-face of a zone.
      REAL xb_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the x-face. Only one is needed
! because this face is flat.
      REAL xb_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the y-face of a zone.
      REAL yb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the y-face of a zone.
      REAL yb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the y-face. Only one is needed
! because this face is flat.
      REAL yb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the z-face of a zone.
      REAL zb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the z-face of a zone.
      REAL zb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the z-face. Because this face is
! curved, each quadrature point will have a unit normal that points in a slightly
! different direction.
      REAL zb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! ----------------------------------------------------------------------
! Variables internal to this subroutine defined here.

      INTEGER ixx, iyy, izz, ixxp1, ixxm1, iyyp1, iyym1, izzm1, i, ii

! "xg, yg, wt, x, y" are needed for 1D quadrature points on the reference element and
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy" are isoparametric mapping functions
! on the reference element. Therefore, we don't need to turn them into vectors.

      REAL x, y, theta, wt1, wt2, wt3, wt4

      REAL, DIMENSION ( NBFLD_QUAD_PTS):: xg, yg

      REAL, DIMENSION ( NBFLD_QUAD_PTS) :: wt

      REAL, DIMENSION ( 8, NBFLD_QUAD_PTS) ::
     1      basis_2D, d_basis_2D_dx, d_basis_2D_dy

! Debug arrays.
      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost ) ::
     1      total_x_area, total_y_area, max_x_normal_deviation,
     1      max_y_normal_deviation

      REAL, DIMENSION ( 0 - NGhost: nr_zones + NGhost) ::
     1      total_z_area, max_z_normal_deviation,
     1      max_z_radius_deviation,
     1      my_radius1, my_radius2_min, my_radius2_max, my_normal

! The remaining temporary variables are needed in the big vectorizing loops.

      REAL radius1, radius2, tempaa, a

      REAL vert ( 8, 3)

      REAL, DIMENSION ( NBFLD_QUAD_PTS) ::
     1      Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag

! ----------------------------------------------------------------------
! STEP 1.a) Build coordinates and weights of quadrature points in the reference
! element. At all orders, up to 4th order, we can get by with 5 points. This is
! the minimum number needed at 3rd and 4th order.

! ------------------------------

      xg(1) = 0.0
      yg(1) = 0.0
      wt(1) = 1.0 / 3.0

      xg(2) = 0.5
      yg(2) = 0.0
      wt(2) = 1.0 / 6.0

      xg(3) = - 0.5
      yg(3) = 0.0
      wt(3) = 1.0 / 6.0

      xg(4) = 0.0
      yg(4) = 0.5
      wt(4) = 1.0 / 6.0

      xg(5) = 0.0
      yg(5) = - 0.5
      wt(5) = 1.0 / 6.0

! --------------------------------------------------
! STEP 1.b) Now use the above to build the basis functions and their derivatives
! in the 2D reference element. We will use isoparametric mapping with second degree
! polynomial nodal bases. This works at all reconstruction orders.

! NBFLD_QUAD_PTS == "5", for 2nd, 3rd and 4th order.

      DO i = 1, NBFLD_QUAD_PTS

        x = xg ( i)
        y = yg ( i)

! ----------

        basis_2D ( 1, i) = (0.5 - x) * (0.5 - y) * (-2.0*x - 2.0*y -1.0)

        d_basis_2D_dx ( 1, i) = - (2.0 * x + y) * (-1.0 + 2.0 * y)

        d_basis_2D_dy ( 1, i) = - (-1.0 + 2.0 * x) * (x + 2.0 * y)

! ----------

        basis_2D ( 2, i) = (0.5 + x) * (0.5 - y) * ( 2.0*x - 2.0*y -1.0)

        d_basis_2D_dx ( 2, i) = - (2.0*x - y) * (-1.0 + 2.0*y)

        d_basis_2D_dy ( 2, i) = - (1.0 + 2.0*x) * (x - 2.0*y)

! ----------

        basis_2D ( 3, i) = (0.5 + x) * (0.5 + y) * ( 2.0*x + 2.0*y -1.0)

        d_basis_2D_dx ( 3, i) = (2.0*x + y) * (1.0 + 2.0*y)

        d_basis_2D_dy ( 3, i) = (1.0 + 2.0*x) * (x + 2.0*y)

! ----------

        basis_2D ( 4, i) = (0.5 - x) * (0.5 + y) * (-2.0*x + 2.0*y -1.0)

        d_basis_2D_dx ( 4, i) = (2.0*x - y) * (1.0 + 2.0*y)

        d_basis_2D_dy ( 4, i) = (-1.0 + 2.0*x) * (x - 2.0*y)

! ----------

        basis_2D ( 5, i) = (1.0 - 4.0*x**2) * (0.5 - y)

        d_basis_2D_dx ( 5, i) = - 8.0 * x * (0.5 - y)

        d_basis_2D_dy ( 5, i) = - 1.0 + 4.0*x**2

! ----------

        basis_2D ( 6, i) = (0.5 + x) * (1.0 - 4.0*y**2)

        d_basis_2D_dx ( 6, i) = 1.0 - 4.0*y**2

        d_basis_2D_dy ( 6, i) = - 8.0 * (0.5 + x) * y

! ----------

        basis_2D ( 7, i) = (1.0 - 4.0*x**2) * (0.5 + y)

        d_basis_2D_dx ( 7, i) = - 8.0 * x * (0.5 + y)

        d_basis_2D_dy ( 7, i) = 1.0 - 4.0*x**2

! ----------

        basis_2D ( 8, i) = (0.5 - x) * (1.0 - 4.0*y**2)

        d_basis_2D_dx ( 8, i) = -1.0 + 4.0*y**2

        d_basis_2D_dy ( 8, i) = - 8.0 * (0.5 - x) * y

! ----------

      END DO   ! End of "DO i = 1, NBFLD_QUAD_PTS".

! ----------
! Eventually the initialization of these small arrays is best done
! on GPU. That is not yet the case here. So we have done it on CPU and
! then copied it into GPU.

!$acc data copyin( basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt )

! ----------------------------------------------------------------------
! STEP 1.c) Initialize un-physical zones and faces with zero or smallnum
! as needed, along x, y and z - faces and then for the volume

! ----------------------------------------
! initializing un-physical zones along x-faces.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, ixxp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch,                                           &
!$OMP&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    imask_active_patch,                                           &
!$acc&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch)        &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, ixxp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------

          ii = MAX0 ( ixx, 1 - NGhost)
          ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)

          IF ( ( imask_active_patch ( ii, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixxp1, iyy, izz) .EQ. 0) ) THEN

            xb_quadpts_patch  ( ixx, iyy, izz, :, :) = 0.0
            xb_quadwt_patch ( ixx, iyy, izz, :)    = 0.0

            xb_quadnorml_patch ( ixx, iyy, izz, :)   = 0.0

          END IF

! ----------

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------
! initializing un-physical zones along y-faces.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, iyyp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch,                                           &
!$OMP&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    imask_active_patch,                                           &
!$acc&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch)        &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, iyyp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------
          ii = MAX0 ( iyy, 1 - NGhost)
          iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)

          IF ( ( imask_active_patch ( ixx, ii, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx, iyyp1, izz) .EQ. 0) ) THEN

            yb_quadpts_patch  ( ixx, iyy, izz, :, :) = 0.0
            yb_quadwt_patch ( ixx, iyy, izz, :)    = 0.0

            yb_quadnorml_patch ( ixx, iyy, izz, :)   = 0.0

          END IF

! ----------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------
! initializing un-physical zones along z-faces.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch,                                           &
!$OMP&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    imask_active_patch,                                           &
!$acc&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)        &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------
! zones are masked along the "x" and "y" directions only and this is
! independent of the radial direction - "izz".

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 0) THEN

            zb_quadpts_patch  ( ixx, iyy, izz, :, :) = 0.0
            zb_quadwt_patch ( ixx, iyy, izz, :)    = 0.0

            zb_quadnorml_patch ( ixx, iyy, izz, :, :)   = 0.0

          END IF

! ----------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 2.a) Using "vert_coords_patch" as input, along with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy" as input, we build
! "xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch" at the flat x-faces.
! This involves a large number of small matrix-multiplies, which can be
! done very efficiently on a GPU.

! By the end of this loop nest, if all the quadrature points within the
! x-face of "xb_quadwt_patch" are added up, we should get a good approximation
! for the area of that x-face.

! The temporaries "radius1, radius2, tempaa, vert,
! Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag" have been turned into vectors
! with "ixx" as the first vector index.


!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    radius1, radius2, tempaa, vert,                               &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)                     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt,                   &
!$OMP&    vert_coords_patch,                                            &
!$OMP&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt,                   &
!$acc&    vert_coords_patch,                                            &
!$acc&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch)

!$acc loop gang vector collapse(3) independent
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    radius1, radius2, tempaa, vert,                               &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

          iyym1 = iyy - 1
          izzm1 = izz - 1

! ------------------------------
! Extract the coordinates for the four vertices for that face.
! Vertices 1, 2, 3, 4 go counter-clockwise in the yz-plane.
! Because the x-face corresponds to a slice of an annulus, we also obtain
! the inner and outer radial distances.

          vert ( 1, :) = vert_coords_patch ( ixx, iyym1, izzm1, :)
          vert ( 2, :) = vert_coords_patch ( ixx, iyy, izzm1, :)
          vert ( 3, :) = vert_coords_patch ( ixx, iyy, izz, :)
          vert ( 4, :) = vert_coords_patch ( ixx, iyym1, izz, :)

          radius1 = SQRT ( vert ( 1, 1)**2 ! inner radius of annulus.
     1                   + vert ( 1, 2)**2
     1                   + vert ( 1, 3)**2 )

          radius2 = SQRT ( vert ( 3, 1)**2  ! outer radius of annulus.
     1                   + vert ( 3, 2)**2
     1                   + vert ( 3, 3)**2 )

! ------------------------------
! Build the midpoints of the four arcs that bound this face. Normalize
! their length to "radius1" or "radius2". Vertex 5 is between vertices 1 and 2.
! Vertex 6 is between vertices 2 and 3. Vertex 7 is between 3 and 4.
! Vertex 8 is between 4 and 1.

          vert ( 5, :) = 0.5 * ( vert ( 1, :) + vert ( 2, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 5, 1)**2
     1                + vert ( 5, 2)**2 + vert ( 5, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert ( 5, :) = vert ( 5, :) * tempaa

          vert ( 6, :) = 0.5 * ( vert ( 2, :) + vert ( 3, :))

          vert ( 7, :) = 0.5 * ( vert ( 3, :) + vert ( 4, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 7, 1)**2
     1                + vert ( 7, 2)**2 + vert ( 7, 3)**2 ) )
          tempaa = radius2 / tempaa
          vert ( 7, :) = vert ( 7, :) * tempaa

          vert ( 8, :) = 0.5 * ( vert ( 4, :) + vert ( 1, :))

! ------------------------------
! Now we extract x,y,z-components from "vert" and do the MATMUL with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".

! We are already in a position to load in the three components of the
! nodal positions in this face, "xb_quadpts_patch ( ixx, iyy, izz, 1:3, :)".

! We also build the vector ( Ax, Ay, Az) which is the partial derivative
! of the position vector w.r.t. the x-coordinate in the reference element.

! We also build the vector ( Bx, By, Bz) which is the partial derivative
! of the position vector w.r.t. the y-coordinate in the reference element.

! ----------
! Extract all the global x-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 1)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "xb_quadpts_patch ( ixx, iyy, izz, 1, :)".

          xb_quadpts_patch ( ixx, iyy, izz, 1, :) =
     1       vert ( 1, 1) * basis_2D ( 1, :)
     1     + vert ( 2, 1) * basis_2D ( 2, :)
     1     + vert ( 3, 1) * basis_2D ( 3, :)
     1     + vert ( 4, 1) * basis_2D ( 4, :)
     1     + vert ( 5, 1) * basis_2D ( 5, :)
     1     + vert ( 6, 1) * basis_2D ( 6, :)
     1     + vert ( 7, 1) * basis_2D ( 7, :)
     1     + vert ( 8, 1) * basis_2D ( 8, :)

          Ax ( :) =
     1       vert ( 1, 1) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dx ( 8, :)

          Bx ( :) =
     1       vert ( 1, 1) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global y-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 2)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "xb_quadpts_patch ( ixx, iyy, izz, 2, :)".

          xb_quadpts_patch ( ixx, iyy, izz, 2, :) =
     1       vert ( 1, 2) * basis_2D ( 1, :)
     1     + vert ( 2, 2) * basis_2D ( 2, :)
     1     + vert ( 3, 2) * basis_2D ( 3, :)
     1     + vert ( 4, 2) * basis_2D ( 4, :)
     1     + vert ( 5, 2) * basis_2D ( 5, :)
     1     + vert ( 6, 2) * basis_2D ( 6, :)
     1     + vert ( 7, 2) * basis_2D ( 7, :)
     1     + vert ( 8, 2) * basis_2D ( 8, :)

          Ay ( :) =
     1       vert ( 1, 2) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dx ( 8, :)

          By ( :) =
     1       vert ( 1, 2) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global z-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 3)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "xb_quadpts_patch ( ixx, iyy, izz, 3, :)".

          xb_quadpts_patch ( ixx, iyy, izz, 3, :) =
     1       vert ( 1, 3) * basis_2D ( 1, :)
     1     + vert ( 2, 3) * basis_2D ( 2, :)
     1     + vert ( 3, 3) * basis_2D ( 3, :)
     1     + vert ( 4, 3) * basis_2D ( 4, :)
     1     + vert ( 5, 3) * basis_2D ( 5, :)
     1     + vert ( 6, 3) * basis_2D ( 6, :)
     1     + vert ( 7, 3) * basis_2D ( 7, :)
     1     + vert ( 8, 3) * basis_2D ( 8, :)

          Az ( :) =
     1       vert ( 1, 3) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dx ( 8, :)

          Bz ( :) =
     1       vert ( 1, 3) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dy ( 8, :)

! ------------------------------
! Now we obtain the cross product of vectors "A" and "B" and store in "C". We
! also obtain its magnitude and store it in "Cmag". Along with the weight "wt",
! that will give us the nodal weights "xb_quadwt_patch". We will also normalize
! the vector "C" so as to obtain the single normal for all the quadrature
! quadrature points in the x-face "xb_quadnorml_patch".

          Cx ( :) = Ay ( :) * Bz ( :) - Az ( :) * By ( :)
          Cy ( :) = Az ( :) * Bx ( :) - Ax ( :) * Bz ( :)
          Cz ( :) = Ax ( :) * By ( :) - Ay ( :) * Bx ( :)

          Cmag ( :) = SQRT ( Cx ( :)**2 + Cy ( :)**2 + Cz ( :)**2)

          xb_quadwt_patch ( ixx, iyy, izz, :) = wt ( :) * Cmag ( :)

          Cmag ( 1) = AMAX1 ( Cmag ( 1), smallnum)

          Cmag ( 1) = 1.0 / Cmag ( 1)

          xb_quadnorml_patch ( ixx, iyy, izz, 1) = Cmag ( 1) * Cx ( 1)
          xb_quadnorml_patch ( ixx, iyy, izz, 2) = Cmag ( 1) * Cy ( 1)
          xb_quadnorml_patch ( ixx, iyy, izz, 3) = Cmag ( 1) * Cz ( 1)

! ------------------------------

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! STEP 2.b) Normalize "xb_quadwt_patch" weights to add up to unity.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, a)                           &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&       xb_quadwt_patch)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present ( smallnum,    &
!$acc&       xb_quadwt_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, a)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

! ------------------------------

          a = xb_quadwt_patch ( ixx, iyy, izz, 1)
     1      + xb_quadwt_patch ( ixx, iyy, izz, 2)
     1      + xb_quadwt_patch ( ixx, iyy, izz, 3)
     1      + xb_quadwt_patch ( ixx, iyy, izz, 4)
     1      + xb_quadwt_patch ( ixx, iyy, izz, 5)

          a = 1.0 / AMAX1 ( a, smallnum)

          xb_quadwt_patch ( ixx, iyy, izz, :) =
     1  a * xb_quadwt_patch ( ixx, iyy, izz, :)

! ------------------------------

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! Make several debugging diagnostics, if that is called for.

#if (DEBUG_GEOM == 1)

! ------------------------------
!$acc update host (                                                     &
!$acc&           vert_coords_patch,                                     &
!$acc&           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch, &
!$acc&           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch, &
!$acc&           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ------------------------------

! We measure how much the normal vector deviates from the true transverse
! direction based on a cross product.
      max_x_normal_deviation = 0.0

      DO iyy = 1 - NGhost, ntp_zones + NGhost
      iyym1 = iyy - 1
        DO izz = 1 - NGhost, nr_zones + NGhost
          DO ixx = 0 - NGhost, ntp_zones + NGhost

          Ax ( 1) = vert_coords_patch ( ixx, iyy, izz, 1)
          Ay ( 1) = vert_coords_patch ( ixx, iyy, izz, 2)
          Az ( 1) = vert_coords_patch ( ixx, iyy, izz, 3)

          Bx ( 1) = vert_coords_patch ( ixx, iyym1, izz, 1)
          By ( 1) = vert_coords_patch ( ixx, iyym1, izz, 2)
          Bz ( 1) = vert_coords_patch ( ixx, iyym1, izz, 3)

          Cx ( 1) = Ay ( 1) * Bz ( 1)
     1                 - Az ( 1) * By ( 1)  ! Cross product of above 2vectors.
          Cy ( 1) = Az ( 1) * Bx ( 1)
     1                 - Ax ( 1) * Bz ( 1)
          Cz ( 1) = Ax ( 1) * By ( 1)
     1                 - Ay ( 1) * Bx ( 1)

          Cmag ( 1) = 1.0 / AMAX1 ( smallnum, SQRT( Cx ( 1)**2
     1                            + Cy ( 1)**2 + Cz ( 1)**2) )

          Cx ( 1) = Cx ( 1) * Cmag ( 1)  ! Normalize the cross product.
          Cy ( 1) = Cy ( 1) * Cmag ( 1)
          Cz ( 1) = Cz ( 1) * Cmag ( 1)

          Cx ( 1) = Cx ( 1)
     1                 - xb_quadnorml_patch ( ixx, iyy, izz, 1) !Diff wrt normal vector.
          Cy ( 1) = Cy ( 1)
     1                 - xb_quadnorml_patch ( ixx, iyy, izz, 2)
          Cz ( 1) = Cz ( 1)
     1                 - xb_quadnorml_patch ( ixx, iyy, izz, 3)

          max_x_normal_deviation ( ixx) =
     1    AMAX1 ( max_x_normal_deviation ( ixx),
     1    SQRT ( Cx ( 1)**2 + Cy ( 1)**2 + Cz ( 1)**2 ) )

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.
        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.
      END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      DO ixx = 0 - NGhost, ntp_zones + NGhost
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ixx, max_x_normal_deviation = ",
     1              ixx, max_x_normal_deviation ( ixx)
      END DO
      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"

! ------------------------------

#endif
! End of DEBUG section.

! ----------------------------------------------------------------------
! STEP 3.a) Using "vert_coords_patch" as input, along with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy" as input, we build
! "yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch" at the flat y-faces.
! This involves a large number of small matrix-multiplies, which can be
! done very efficiently on a GPU.

! By the end of this loop nest, if all the quadrature points within the
! y-face of "yb_quadwt_patch" are added up, we should get a good approximation
! for the area of that y-face.

! The temporaries "radius1, radius2, tempaa, vert,
! Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag" have been turned into vectors
! with "iyy" as the first vector index.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    radius1, radius2, tempaa, vert,                               &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)                     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt,                   &

!$OMP&    vert_coords_patch,                                            &
!$OMP&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt,                   &
!$acc&    vert_coords_patch,                                            &
!$acc&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    radius1, radius2, tempaa, vert,                               &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

          izzm1 = izz - 1
          ixxm1 = ixx - 1

! ------------------------------
! Extract the coordinates for the four vertices for that face.
! Vertices 1, 2, 3, 4 go counter-clockwise in the zx-plane.
! Because the y-face corresponds to a slice of an annulus, we also obtain
! the inner and outer radial distances.

          vert ( 1, :) = vert_coords_patch ( ixxm1, iyy, izzm1, :)
          vert ( 2, :) = vert_coords_patch ( ixxm1, iyy, izz, :)
          vert ( 3, :) = vert_coords_patch ( ixx, iyy, izz, :)
          vert ( 4, :) = vert_coords_patch ( ixx, iyy, izzm1, :)

          radius1 = SQRT ( vert ( 1, 1)**2 ! inner radius of annulus.
     1                   + vert ( 1, 2)**2
     1                   + vert ( 1, 3)**2 )

          radius2 = SQRT ( vert ( 3, 1)**2  ! outer radius of annulus.
     1                   + vert ( 3, 2)**2
     1                   + vert ( 3, 3)**2 )

! ------------------------------
! Build the midpoints of the four arcs that bound this face. Normalize
! their length to "radius1" or "radius2". Vertex 5 is between vertices 1 and 2.
! Vertex 6 is between vertices 2 and 3. Vertex 7 is between 3 and 4.
! Vertex 8 is between 4 and 1.

          vert( 5, :) = 0.5 * ( vert( 1, :) + vert( 2, :))

          vert( 6, :) = 0.5 * ( vert( 2, :) + vert( 3, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert( 6, 1)**2
     1                   + vert( 6, 2)**2 + vert( 6, 3)**2 ) )
          tempaa = radius2 / tempaa
          vert( 6, :) = vert( 6, :) * tempaa

          vert( 7, :) = 0.5 * ( vert( 3, :) + vert( 4, :))

          vert( 8, :) = 0.5 * ( vert( 4, :) + vert( 1, :))

          tempaa = AMAX1 ( smallnum, SQRT ( vert( 8, 1)**2
     1                   + vert( 8, 2)**2 + vert( 8, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert( 8, :) = vert( 8, :) * tempaa

! ------------------------------
! Now we extract x,y,z-components from "vert" and do the MATMUL with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".

! We are already in a position to load in the three components of the
! nodal positions in this face, "yb_quadpts_patch ( ixx, iyy, izz, 1:3, :)".

! We also build the vector ( Ax, Ay, Az) which is the partial derivative
! of the position vector w.r.t. the x-coordinate in the reference element.

! We also build the vector ( Bx, By, Bz) which is the partial derivative
! of the position vector w.r.t. the y-coordinate in the reference element.

! ----------
! Extract all the global x-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 1)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "yb_quadpts_patch ( ixx, iyy, izz, 1, :)".

          yb_quadpts_patch ( ixx, iyy, izz, 1, :) =
     1       vert ( 1, 1) * basis_2D ( 1, :)
     1     + vert ( 2, 1) * basis_2D ( 2, :)
     1     + vert ( 3, 1) * basis_2D ( 3, :)
     1     + vert ( 4, 1) * basis_2D ( 4, :)
     1     + vert ( 5, 1) * basis_2D ( 5, :)
     1     + vert ( 6, 1) * basis_2D ( 6, :)
     1     + vert ( 7, 1) * basis_2D ( 7, :)
     1     + vert ( 8, 1) * basis_2D ( 8, :)

          Ax ( :) =
     1       vert ( 1, 1) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dx ( 8, :)

          Bx ( :) =
     1       vert ( 1, 1) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global y-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 2)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "xb_quadpts_patch ( ixx, iyy, izz, 2, :)".

          yb_quadpts_patch ( ixx, iyy, izz, 2, :) =
     1       vert ( 1, 2) * basis_2D ( 1, :)
     1     + vert ( 2, 2) * basis_2D ( 2, :)
     1     + vert ( 3, 2) * basis_2D ( 3, :)
     1     + vert ( 4, 2) * basis_2D ( 4, :)
     1     + vert ( 5, 2) * basis_2D ( 5, :)
     1     + vert ( 6, 2) * basis_2D ( 6, :)
     1     + vert ( 7, 2) * basis_2D ( 7, :)
     1     + vert ( 8, 2) * basis_2D ( 8, :)

          Ay ( :) =
     1       vert ( 1, 2) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dx ( 8, :)

          By ( :) =
     1       vert ( 1, 2) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global z-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 3)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "yb_quadpts_patch ( ixx, iyy, izz, 3, :)".

          yb_quadpts_patch ( ixx, iyy, izz, 3, :) =
     1       vert ( 1, 3) * basis_2D ( 1, :)
     1     + vert ( 2, 3) * basis_2D ( 2, :)
     1     + vert ( 3, 3) * basis_2D ( 3, :)
     1     + vert ( 4, 3) * basis_2D ( 4, :)
     1     + vert ( 5, 3) * basis_2D ( 5, :)
     1     + vert ( 6, 3) * basis_2D ( 6, :)
     1     + vert ( 7, 3) * basis_2D ( 7, :)
     1     + vert ( 8, 3) * basis_2D ( 8, :)

          Az ( :) =
     1       vert ( 1, 3) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dx ( 8, :)

          Bz ( :) =
     1       vert ( 1, 3) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dy ( 8, :)

! ------------------------------
! Now we obtain the cross product of vectors "A" and "B" and store in "C". We
! also obtain its magnitude and store it in "Cmag". Along with the weight "wt",
! that will give us the nodal weights "yb_quadwt_patch". We will also normalize
! the vector "C" so as to obtain directions of the normals at each quadrature
! point "yb_quadnorml_patch".

          Cx ( :) = Ay ( :) * Bz ( :) - Az ( :) * By ( :)
          Cy ( :) = Az ( :) * Bx ( :) - Ax ( :) * Bz ( :)
          Cz ( :) = Ax ( :) * By ( :) - Ay ( :) * Bx ( :)

          Cmag ( :) = SQRT ( Cx ( :)**2 + Cy ( :)**2 + Cz ( :)**2)

          yb_quadwt_patch ( ixx, iyy, izz, :) = wt ( :) * Cmag ( :)

          Cmag ( 1) = AMAX1 ( Cmag ( 1), smallnum)

          Cmag ( 1) = 1.0 / Cmag ( 1)

          yb_quadnorml_patch ( ixx, iyy, izz, 1) = Cmag ( 1) * Cx ( 1)
          yb_quadnorml_patch ( ixx, iyy, izz, 2) = Cmag ( 1) * Cy ( 1)
          yb_quadnorml_patch ( ixx, iyy, izz, 3) = Cmag ( 1) * Cz ( 1)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! STEP 3.b) Normalize "yb_quadwt_patch" weights to add up to unity.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, a)                           &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    yb_quadwt_patch)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present ( smallnum,    &
!$acc&    yb_quadwt_patch)                                              &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, a)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

! ------------------------------

          a = yb_quadwt_patch ( ixx, iyy, izz, 1)
     1      + yb_quadwt_patch ( ixx, iyy, izz, 2)
     1      + yb_quadwt_patch ( ixx, iyy, izz, 3)
     1      + yb_quadwt_patch ( ixx, iyy, izz, 4)
     1      + yb_quadwt_patch ( ixx, iyy, izz, 5)

          a = 1.0 / AMAX1 ( a, smallnum)

          yb_quadwt_patch ( ixx, iyy, izz, :) =
     1    a * yb_quadwt_patch ( ixx, iyy, izz, :)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

!--------------------------------------------------
! Make several debugging diagnostics, if that is called for.

#if (DEBUG_GEOM == 1)

! ------------------------------
!$acc update host (                                                     &
!$acc&           vert_coords_patch,                                     &
!$acc&           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch, &
!$acc&           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch, &
!$acc&           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ------------------------------

! We measure how much the normal vector deviates from the true transverse
! direction based on a cross product.
      max_y_normal_deviation = 0.0

      DO ixx = 1 - NGhost, ntp_zones + NGhost
      ixxm1 = ixx - 1
        DO izz = 1 - NGhost, nr_zones + NGhost
          DO iyy = 0 - NGhost, ntp_zones + NGhost

          Ax ( 1) = vert_coords_patch ( ixxm1, iyy, izz, 1)
          Ay ( 1) = vert_coords_patch ( ixxm1, iyy, izz, 2)
          Az ( 1) = vert_coords_patch ( ixxm1, iyy, izz, 3)

          Bx ( 1) = vert_coords_patch ( ixx, iyy, izz, 1)
          By ( 1) = vert_coords_patch ( ixx, iyy, izz, 2)
          Bz ( 1) = vert_coords_patch ( ixx, iyy, izz, 3)

          Cx ( 1) = Ay ( 1) * Bz ( 1)
     1                 - Az ( 1) * By ( 1)  ! Cross product of above 2 vectors.
          Cy ( 1) = Az ( 1) * Bx ( 1)
     1                 - Ax ( 1) * Bz ( 1)
          Cz ( 1) = Ax ( 1) * By ( 1)
     1                 - Ay ( 1) * Bx ( 1)

          Cmag ( 1) = 1.0 / AMAX1 ( smallnum, SQRT( Cx ( 1)**2
     1              + Cy ( 1)**2 + Cz ( 1)**2) )

          Cx ( 1) = Cx ( 1) * Cmag ( 1)  ! Normalize the cross product.
          Cy ( 1) = Cy ( 1) * Cmag ( 1)
          Cz ( 1) = Cz ( 1) * Cmag ( 1)

          Cx ( 1) = Cx ( 1)
     1                 - yb_quadnorml_patch ( ixx, iyy, izz, 1)  ! Diff wrt normal vector.
          Cy ( 1) = Cy ( 1)
     1                 - yb_quadnorml_patch ( ixx, iyy, izz, 2)
          Cz ( 1) = Cz ( 1)
     1                 - yb_quadnorml_patch ( ixx, iyy, izz, 3)

          max_y_normal_deviation ( iyy) =
     1          AMAX1 ( max_y_normal_deviation ( iyy),
     1    SQRT ( Cx ( 1)**2 + Cy ( 1)**2 + Cz ( 1)**2 ) )

          END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.
        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.
      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

      DO iyy = 0 - NGhost, ntp_zones + NGhost
      WRITE ( 6, *)" "
      WRITE ( 6, *)"iyy, max_y_normal_deviation = ",
     1              iyy, max_y_normal_deviation ( iyy)
      END DO
      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"

! ------------------------------

#endif
! End of DEBUG part.

! ----------------------------------------------------------------------
! STEP 4.a) Using "vert_coords_patch" as input, along with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy" as input, we build
! "zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch" at the curved z-faces.
! This involves a large number of small matrix-multiplies, which can be
! done very efficiently on a GPU.

! By the end of this loop nest, if all the quadrature points within the
! z-face of "zb_quadwt_patch" are added up, we should get a good approximation
! for the area of that z-face.

! The temporaries "radius1, radius2, tempaa, vert,
! Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag" have been turned into vectors
! with "izz" as the first vector index.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    radius1, radius2, tempaa, vert,                               &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)                     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt, smallnum,         &

!$OMP&    vert_coords_patch,                                            &
!$OMP&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt,                   &
!$acc&    vert_coords_patch,                                            &
!$acc&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    radius1, radius2, tempaa, vert,                               &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Cmag)


      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------
! Extract the coordinates for the four vertices for that face.
! Vertices 1, 2, 3, 4 go counter-clockwise in the xy-plane.
! Because the z-direction is the radial direction, we also obtain the radial distance.

          vert ( 1, :) = vert_coords_patch ( ixxm1, iyym1, izz, :)
          vert ( 2, :) = vert_coords_patch ( ixx, iyym1, izz, :)
          vert ( 3, :) = vert_coords_patch ( ixx, iyy, izz, :)
          vert ( 4, :) = vert_coords_patch ( ixxm1, iyy, izz, :)

          radius1 = SQRT ( vert ( 1, 1)**2
     1                   + vert ( 1, 2)**2  ! Radius of spherical patch.
     1                   + vert ( 1, 3)**2 )

! ------------------------------
! Build the midpoints of the four arcs that bound this face. Normalize
! their length to "radius1". Vertex 5 is between vertices 1 and 2.
! Vertex 6 is between vertices 2 and 3. Vertex 7 is between 3 and 4.
! Vertex 8 is between 4 and 1.

          vert ( 5, :) = 0.5 * ( vert ( 1, :) + vert ( 2, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 5, 1)**2
     1                 + vert ( 5, 2)**2 + vert ( 5, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert ( 5, :) = vert ( 5, :) * tempaa

! ----------
          vert ( 6, :) = 0.5 * ( vert ( 2, :) + vert ( 3, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 6, 1)**2
     1                 + vert ( 6, 2)**2 + vert ( 6, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert ( 6, :) = vert ( 6, :) * tempaa

! ----------
          vert ( 7, :) = 0.5 * ( vert ( 3, :) + vert ( 4, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 7, 1)**2
     1                 + vert ( 7, 2)**2 + vert ( 7, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert ( 7, :) = vert ( 7, :) * tempaa

! ----------
          vert ( 8, :) = 0.5 * ( vert ( 4, :) + vert ( 1, :))
          tempaa = AMAX1 ( smallnum, SQRT ( vert ( 8, 1)**2
     1                 + vert ( 8, 2)**2 + vert ( 8, 3)**2 ) )
          tempaa = radius1 / tempaa
          vert ( 8, :) = vert ( 8, :) * tempaa

! ------------------------------
! Now we extract x,y,z-components from "vert" and do the MATMUL with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".

! We are already in a position to load in the three components of the
! nodal positions in this face, "zb_quadpts_patch ( ixx, iyy, izz, 1:3, :)".

! We also build the vector ( Ax, Ay, Az) which is the partial derivative
! of the position vector w.r.t. the x-coordinate in the reference element.

! We also build the vector ( Bx, By, Bz) which is the partial derivative
! of the position vector w.r.t. the y-coordinate in the reference element.

! ----------
! Extract all the global x-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 1)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "zb_quadpts_patch ( ixx, iyy, izz, 1, :)".

          zb_quadpts_patch ( ixx, iyy, izz, 1, :) =
     1       vert ( 1, 1) * basis_2D ( 1, :)
     1     + vert ( 2, 1) * basis_2D ( 2, :)
     1     + vert ( 3, 1) * basis_2D ( 3, :)
     1     + vert ( 4, 1) * basis_2D ( 4, :)
     1     + vert ( 5, 1) * basis_2D ( 5, :)
     1     + vert ( 6, 1) * basis_2D ( 6, :)
     1     + vert ( 7, 1) * basis_2D ( 7, :)
     1     + vert ( 8, 1) * basis_2D ( 8, :)

          Ax ( :) =
     1       vert ( 1, 1) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dx ( 8, :)

          Bx ( :) =
     1       vert ( 1, 1) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 1) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 1) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 1) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 1) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 1) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 1) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 1) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global y-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 2)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "xb_quadpts_patch ( ixx, iyy, izz, 2, :)".

          zb_quadpts_patch ( ixx, iyy, izz, 2, :) =
     1       vert ( 1, 2) * basis_2D ( 1, :)
     1     + vert ( 2, 2) * basis_2D ( 2, :)
     1     + vert ( 3, 2) * basis_2D ( 3, :)
     1     + vert ( 4, 2) * basis_2D ( 4, :)
     1     + vert ( 5, 2) * basis_2D ( 5, :)
     1     + vert ( 6, 2) * basis_2D ( 6, :)
     1     + vert ( 7, 2) * basis_2D ( 7, :)
     1     + vert ( 8, 2) * basis_2D ( 8, :)

          Ay ( :) =
     1       vert ( 1, 2) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dx ( 8, :)

          By ( :) =
     1       vert ( 1, 2) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 2) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 2) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 2) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 2) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 2) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 2) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 2) * d_basis_2D_dy ( 8, :)

! ----------
! Extract all the global z-coordinates and their x and y derivatives w.r.t.
! the reference element by multiplying "vert ( 1:8, 3)" with
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy".
! Store in "zb_quadpts_patch ( ixx, iyy, izz, 3, :)".

          zb_quadpts_patch ( ixx, iyy, izz, 3, :) =
     1       vert ( 1, 3) * basis_2D ( 1, :)
     1     + vert ( 2, 3) * basis_2D ( 2, :)
     1     + vert ( 3, 3) * basis_2D ( 3, :)
     1     + vert ( 4, 3) * basis_2D ( 4, :)
     1     + vert ( 5, 3) * basis_2D ( 5, :)
     1     + vert ( 6, 3) * basis_2D ( 6, :)
     1     + vert ( 7, 3) * basis_2D ( 7, :)
     1     + vert ( 8, 3) * basis_2D ( 8, :)

          Az ( :) =
     1       vert ( 1, 3) * d_basis_2D_dx ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dx ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dx ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dx ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dx ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dx ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dx ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dx ( 8, :)

          Bz ( :) =
     1       vert ( 1, 3) * d_basis_2D_dy ( 1, :)
     1     + vert ( 2, 3) * d_basis_2D_dy ( 2, :)
     1     + vert ( 3, 3) * d_basis_2D_dy ( 3, :)
     1     + vert ( 4, 3) * d_basis_2D_dy ( 4, :)
     1     + vert ( 5, 3) * d_basis_2D_dy ( 5, :)
     1     + vert ( 6, 3) * d_basis_2D_dy ( 6, :)
     1     + vert ( 7, 3) * d_basis_2D_dy ( 7, :)
     1     + vert ( 8, 3) * d_basis_2D_dy ( 8, :)

! ------------------------------
! Now we obtain the cross product of vectors "A" and "B" and store in "C". We
! also obtain its magnitude and store it in "Cmag". Along with the weight "wt",
! that will give us the nodal weights "zb_quadwt_patch". We will also normalize
! the vector "C" so as to obtain directions of the normals at each quadrature
! point "zb_quadnorml_patch".

          Cx ( :) = Ay ( :) * Bz ( :) - Az ( :) * By ( :)
          Cy ( :) = Az ( :) * Bx ( :) - Ax ( :) * Bz ( :)
          Cz ( :) = Ax ( :) * By ( :) - Ay ( :) * Bx ( :)

          Cmag ( :) = SQRT ( Cx ( :)**2 + Cy ( :)**2  + Cz ( :)**2)

          zb_quadwt_patch ( ixx, iyy, izz, :) = wt ( :) * Cmag ( :)

          Cmag ( :) = AMAX1 ( Cmag ( :), smallnum)

          Cmag ( :) = 1.0 / Cmag ( :)

          zb_quadnorml_patch ( ixx, iyy, izz, 1, :) =
     1                            Cmag ( :) * Cx ( :)
          zb_quadnorml_patch ( ixx, iyy, izz, 2, :) =
     1                            Cmag ( :) * Cy ( :)
          zb_quadnorml_patch ( ixx, iyy, izz, 3, :) =
     1                            Cmag ( :) * Cz ( :)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! STEP 4.b) Normalize "zb_quadwt_patch" weights to add up to unity.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, a)                           &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    zb_quadwt_patch)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present ( smallnum,    &
!$acc&    zb_quadwt_patch)                                              &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, a)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost  ! This loop can be vectorized.

! ------------------------------

          a = zb_quadwt_patch ( ixx, iyy, izz, 1)
     1      + zb_quadwt_patch ( ixx, iyy, izz, 2)
     1      + zb_quadwt_patch ( ixx, iyy, izz, 3)
     1      + zb_quadwt_patch ( ixx, iyy, izz, 4)
     1      + zb_quadwt_patch ( ixx, iyy, izz, 5)

          a = 1.0 / AMAX1 ( a, smallnum)

          zb_quadwt_patch ( ixx, iyy, izz, :) =
     1  a * zb_quadwt_patch ( ixx, iyy, izz, :)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! --------------------------------------------------
! Make several debugging diagnostics, if that is called for.

#if (DEBUG_GEOM == 1)

! ------------------------------
!$acc update host (                                                     &
!$acc&           vert_coords_patch,                                     &
!$acc&           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch, &
!$acc&           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch, &
!$acc&           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ------------------------------

! We measure how much the normal vector deviates from the true radial direction.
      max_z_normal_deviation = 0.0

! We also measure how much the radius at each of the nodes deviates from the
! radius of a sphere.
      max_z_radius_deviation = 0.0

      my_radius2_min = 1.0e32
      my_radius2_max = -1.0e32

      DO i = 1, NBFLD_QUAD_PTS

      DO ixx = 1 - NGhost, ntp_zones + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
          DO izz = 0 - NGhost, nr_zones + NGhost

! ----------

          Ax ( 1:3) = zb_quadpts_patch ( ixx, iyy, izz, 1:3, i)  ! Radius from isoparm mapping.
          vert ( 1, :) = vert_coords_patch ( ixx, iyy, izz, :)  ! This contains true radius.

          radius1 = SQRT ( vert ( 1, 1)**2
     1                   + vert ( 1, 2)**2  ! Evaluate true radius.
     1                   + vert ( 1, 3)**2 )

          radius2 = SQRT ( Ax ( 1)**2 + Ax ( 2)**2  ! Evaluate radius from isoparm mapping.
     1                   + Ax ( 3)**2 )

          my_radius1 ( izz) = radius1
          my_radius2_min ( izz)= AMIN1(my_radius2_min( izz), radius2)
          my_radius2_max ( izz)= AMAX1(my_radius2_max( izz), radius2)

          max_z_radius_deviation ( izz) =
     1                 AMAX1 ( max_z_radius_deviation ( izz),  ! Deviation in radii.
     1                             ABS ( radius1 - radius2 ) )

! ----------

          Ax ( 1:3) = Ax ( 1:3) / radius2 ! This is truly a unit vector in radial direction.
          Bx ( 1:3) = zb_quadnorml_patch ( ixx, iyy, izz, 1:3, i)  ! Normal from isoparm map.

          my_normal = SQRT ( Bx ( 1)**2 + Bx ( 2)**2
     1                     + Bx ( 3)**2 )

          max_z_normal_deviation ( izz) =
     1          AMAX1 ( max_z_normal_deviation ( izz),  ! Deviation in unit normal.
     1 SQRT ( ( Ax ( 1) - Bx ( 1))**2
     1      + ( Ax ( 2) - Bx ( 2))**2
     1      + ( Ax ( 3) - Bx ( 3))**2 ) )

! ----------

          END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.
        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.
      END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO i = 1, NBFLD_QUAD_PTS" loop.


      DO izz = 0 - NGhost, nr_zones + NGhost
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1      "izz, max_z_normal_deviation, max_z_radius_deviation = ",
     1 izz, max_z_normal_deviation ( izz), max_z_radius_deviation ( izz)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"error ratio = ", max_z_radius_deviation ( izz)
     1                       / my_radius1 ( izz)
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1 "my_radius1, my_radius2_min, max, my_normal length = ",
     1  my_radius1 ( izz), my_radius2_min ( izz), my_radius2_max ( izz),
     1  my_normal ( izz)
      END DO
      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"

! ------------------------------

#endif
! End of DEBUG part.

!$acc end data
! This is because we are releasing the allocation on :
! "basis_2D, d_basis_2D_dx, d_basis_2D_dy, wt".

! ----------------------------------------------------------------------

      END SUBROUTINE COORDS_TO_WT_QUADPTS_BFLD_ISO2


