








      SUBROUTINE GENERATE_GLOBAL_RADIAL_UNIFORM_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Generate radial binning. Logarithmic ratioing set by default, but any
! radial meshing can be initialized here.

! "rmin, rmax" are global min and max values of the radial extent ( input)

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input) on the global AMR mesh.

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone. This is a simple
!  arithmetic averaging of the boundary values. So it is not really a centroid;
!  therefore, do not take it too seriously.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: rmin, rmax, smallnum, radial_mesh_ratio
      INTEGER :: nr_global, NGhost_global

      REAL
     1 rb_save_global ( 0 - NGhost_global: nr_global + NGhost_global)
      REAL
     1 dr_save_global ( 1  - NGhost_global: nr_global + NGhost_global)
      REAL
     1 rc_save_global ( 1  - NGhost_global: nr_global + NGhost_global)

      INTEGER :: i
      REAL :: pow

! ----------------------------------------------------------------------

        DO i = 0 - NGhost_global, nr_global + NGhost_global

          pow = REAL ( i) / REAL ( nr_global)
          rb_save_global ( i) = rmin + ( rmax - rmin)*pow

        END DO

        DO i = 1 - NGhost_global, nr_global + NGhost_global
          dr_save_global ( i) = rb_save_global ( i)
     1                        - rb_save_global ( i - 1)
          rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                                + rb_save_global ( i - 1) )
        END DO

! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_UNIFORM_MESH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GENERATE_GLOBAL_RADIAL_LOG_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Generate radial binning. Logarithmic ratioing set by default, but any
! radial meshing can be initialized here.

! "rmin, rmax" are global min and max values of the radial extent ( input)

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input) on the global AMR mesh.

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone. This is a simple
!  arithmetic averaging of the boundary values. So it is not really a centroid;
!  therefore, do not take it too seriously.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: rmin, rmax, smallnum, radial_mesh_ratio
      INTEGER :: nr_global, NGhost_global

      REAL
     1 rb_save_global ( 0 - NGhost_global: nr_global + NGhost_global)
      REAL
     1 dr_save_global ( 1  - NGhost_global: nr_global + NGhost_global)
      REAL
     1 rc_save_global ( 1  - NGhost_global: nr_global + NGhost_global)

      INTEGER :: i
      REAL :: pow

! ----------------------------------------------------------------------

        DO i = 0 - NGhost_global, nr_global + NGhost_global

          pow = REAL ( i) / REAL ( nr_global)
          rb_save_global ( i) = rmin * ( rmax / rmin)**pow

        END DO

        DO i = 1 - NGhost_global, nr_global + NGhost_global
          dr_save_global ( i) = rb_save_global ( i)
     1                        - rb_save_global ( i - 1)
          rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                                + rb_save_global ( i - 1) )
        END DO

! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_LOG_MESH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GENERATE_GLOBAL_RADIAL_STRETCHED_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Generate geometrically stretched mesh.
!         rmax - rmin = dx0 * ( ratio^N  - 1) / (ratio - 1)

! "rmin, rmax" are global min and max values of the radial extent ( input)
! "radial_mesh_ratio" ratio of the mesh in the radial direction. Set to 1.0 for unratioed.

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input)

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: rmin, rmax, smallnum, radial_mesh_ratio
      INTEGER :: nr_global, NGhost_global

      REAL
     1 rb_save_global ( 0 - NGhost_global: nr_global + NGhost_global)
      REAL
     1 dr_save_global ( 1  - NGhost_global: nr_global + NGhost_global)
      REAL
     1 rc_save_global ( 1  - NGhost_global: nr_global + NGhost_global)

      INTEGER :: i, N
      REAL :: dx0, dx1, dxN

! ----------------------------------------------------------------------

      N = nr_global + 1

! check for meaningful ratio.
      IF ( radial_mesh_ratio .GT. ( 1.0 + smallnum)) THEN

        dx1 = ( rmax - rmin) / ( ( radial_mesh_ratio**N - 1) /
     1                           ( radial_mesh_ratio - 1))

! else uniform mesh.
      ELSE

        dx1 = ( rmax - rmin) / nr_global

        radial_mesh_ratio = 1.0

      END IF

! --------------------

        rb_save_global ( 0) = rmin

        DO i = 1, nr_global

          rb_save_global ( i) = dx1 * radial_mesh_ratio**( i - 1)
     1                              + rb_save_global( i - 1)

        END DO

! ----------

        dx0 = rb_save_global ( 1) - rb_save_global ( 0)

        DO i = 0 - NGhost_global, -1

          rb_save_global ( i) = rb_save_global ( 0) + ( i * dx0)

        END DO

! ----------

        dxN = rb_save_global ( nr_global)
     1      - rb_save_global ( nr_global - 1)

        DO i = nr_global + 1, nr_global + NGhost_global

          rb_save_global ( i) = rb_save_global ( nr_global)
     1                 + ( ( i - nr_global) * dxN)

        END DO

! ----------

        DO i = 1 - NGhost_global, nr_global + NGhost_global

          dr_save_global ( i) = rb_save_global ( i)
     1                        - rb_save_global ( i - 1)

          rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                                + rb_save_global ( i - 1) )

        END DO

! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_STRETCHED_MESH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GENERATE_GLOBAL_RADIAL_HYBRID_MESH

     1       ( nr_zones, nr_global, NGhost_global, nr_geom_zones,
     1         rmin, rmax, smallnum, radial_mesh_ratio, dr1,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Generate hybrid radial mesh.

! "rmin, rmax" are global min and max values of the radial extent ( input)
! "radial_mesh_ratio" ratio of the mesh in the radial direction.
! Set to 1.0 for unratioed.

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input)

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: nr_zones, nr_global, NGhost_global, nr_geom_zones
      REAL :: rmin, rmax, smallnum, radial_mesh_ratio, dr1

      REAL rb_save_global ( 0 - NGhost_global:nr_global + NGhost_global)
      REAL dr_save_global ( 1 - NGhost_global:nr_global + NGhost_global)
      REAL rc_save_global ( 1 - NGhost_global:nr_global + NGhost_global)

! ----------------------------------------------------------------------
! Variables below are local to this subroutine.

      REAL rb_log ( 0 - NGhost_global: nr_global
     1                - nr_geom_zones + NGhost_global)
      REAL dr_log ( 1 - NGhost_global: nr_global
     1                - nr_geom_zones + NGhost_global)
      REAL rc_log ( 1 - NGhost_global: nr_global
     1                - nr_geom_zones + NGhost_global)

! ----------

      REAL rb_geo ( 0 - NGhost_global: nr_geom_zones
     1                + NGhost_global)

      REAL dr_geo ( 1 - NGhost_global: nr_geom_zones
     1                + NGhost_global)

      REAL rc_geo ( 1 - NGhost_global: nr_geom_zones
     1                + NGhost_global)

! ----------

      INTEGER :: i, nr_global_logmesh
      REAL :: pow, drN, dr0

      REAL :: r_geomin, r_geomax, r_logmin, r_logmax, ratio_geom

! ----------------------------------------------------------------------

      ratio_geom = radial_mesh_ratio

! ----------

      r_geomin = rmin

      rb_geo(0) = r_geomin

      DO i = 1, nr_geom_zones

        rb_geo ( i) = dr1 * ratio_geom**( i - 1) + rb_geo( i - 1)

      END DO

! ----------------------
! Compute log mesh on the Radial Slabs 2 to Num_Radial_Slabs

      nr_global_logmesh = nr_global - nr_geom_zones

      r_geomax = rb_geo ( nr_geom_zones)

      r_logmin = r_geomax
      r_logmax = rmax

      DO i = 0, nr_global_logmesh + NGhost_global

        pow = REAL ( i) / REAL ( nr_global_logmesh)
        rb_log ( i) = r_logmin * ( r_logmax / r_logmin)**pow

      END DO

! --------------------
! Filling "rb_save_global" with "rb_geo" and "rb_log".

! filling the inner ghost zones with geometric mesh
      dr0 = rb_geo ( 1) - rb_geo ( 0)

      DO i = 0 - NGhost_global, -1

        rb_save_global ( i) = rb_geo ( 0) + ( i * dr0)

      END DO

! filling 0 - nr_zones with geometric mesh
      DO i = 0, nr_geom_zones

        rb_save_global ( i) = rb_geo ( i)

      END DO

! filling nr_zones + 1 till end with log ratioed mesh
      DO i = nr_geom_zones + 1, nr_global + NGhost_global

        rb_save_global ( i) = rb_log ( i - nr_geom_zones)

      END DO

! Compute "dr_save_global" and "rc_save_global".
      DO i = 1 - NGhost_global, nr_global + NGhost_global

        dr_save_global ( i) = rb_save_global ( i)
     1                      - rb_save_global ( i - 1)

        rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                              + rb_save_global ( i - 1) )

      END DO

! ----------------------------------------------------------------------
! Code to printout and check, do only 1-processor run.
#if 0
        DO i = 1 - NGhost_global, nr_global + NGhost_global

          WRITE( 6, *)"i, rb_save_global ( i), dr_save_global ( i) = ",
     1                 i, rb_save_global ( i), dr_save_global ( i)

        END DO
        CALL FLUSH ( 6)
        ERROR STOP
#endif
! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_HYBRID_MESH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GENERATE_GLOBAL_RADIAL_POWERLAW_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Written by : Vladimir Florinski

! Uses "rmin, rmax" and "nr_global" to provide an increasing mesh spacing
! with increasing radius.
! "radial_mesh_ratio << 1" gives meshes that are smaller at "rmax".
! "radial_mesh_ratio == 1" gives uniform mesh in radial direction.
! "radial_mesh_ratio --> Infinity" gives logarithmic scaling.

! "rmin, rmax" are global min and max values of the radial extent ( input)

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input) on the global AMR mesh.

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone. This is a simple
!  arithmetic averaging of the boundary values. So it is not really a centroid;
!  therefore, do not take it too seriously.

! "radial_mesh_ratio" This is the factor b in equation 26 of our notess.
! When b=1 it gives uniform mesh grid and when b is a very high number,say infinity, it
! reduces to a log scales. For b= 1.0e10 the nonuniform mesh grid points are
! consistent upto the order 10^-4 with the logarithmic mesh grid points.
! Thus increasing the parameter b leads to smaller radial zones close to the origin and
! larger radial zones towards the outer boundary.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: rmin, rmax, smallnum, radial_mesh_ratio
      INTEGER :: nr_global, NGhost_global

      REAL
     1 rb_save_global ( 0 - NGhost_global: nr_global + NGhost_global)
      REAL
     1 dr_save_global ( 1  - NGhost_global: nr_global + NGhost_global)
      REAL
     1 rc_save_global ( 1  - NGhost_global: nr_global + NGhost_global)

      INTEGER :: i
      REAL :: b, dx0, dxN

! ----------------------------------------------------------------------

      b = radial_mesh_ratio

      DO i = 0, nr_global

        rb_save_global ( i) = rmin * ( 1.0
     1          + ( (rmax/rmin)**(1.0/b) - 1.0)
     1          * (REAL (i)/REAL (nr_global))) ** b

      END DO

! ----------

      dx0 = rb_save_global ( 1) - rb_save_global ( 0)

      DO i = 0 - NGhost_global, -1

        rb_save_global ( i) = rb_save_global ( 0) + ( i * dx0)

      END DO

! ----------

      dxN = rb_save_global ( nr_global)
     1    - rb_save_global ( nr_global - 1)

      DO i = nr_global + 1, nr_global + NGhost_global

        rb_save_global ( i) = rb_save_global ( nr_global)
     1               + ( ( i - nr_global) * dxN)

      END DO

! ----------

      DO i = 1 - NGhost_global, nr_global + NGhost_global

        dr_save_global ( i) = rb_save_global ( i)
     1                      - rb_save_global ( i - 1)

        rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                              + rb_save_global ( i - 1) )

      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_POWERLAW_MESH

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GENERATE_GLOBAL_RADIAL_UNIF_N_LOG_MESH
     1       ( nr_global, NGhost_global, nr_geom_zones,
     1         rmin, rmax, rint, smallnum,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Generate radial mesh. It will have total zones "nr_global" out of which
! the inner "nr_global - nr_geom_zones" will be uniform and the outer
! "nr_geom_zones" will be log ratioed. "rint" is the location where this
! switch is made.

! Code does not ensure that the zones will be smoothly varying at "rint".
! Therefore this has to be done in a separate calculation.

! The mesh function for the inner zones is
! "f (xi) = rmin + xi * ((rint - rmin)/(nr_global - nr_geom_zones))"
! with "xi" going from [0, nr_global - nr_geom_zones].

! The mesh function for the outer log ratioed zones is "g (psi) = rint * EXP ( alpha * psi)"
! with "psi" going from [0, nr_geom_zones].

! For this to work, we set
! "alpha = LN ( 1.0 + (rint - rmin)/(rint * (nr_global - nr_geom_zones)) )".

! "rmax" will then have to be reset as "rint * EXP (alpha * nr_geom_zones)".

! "rmin, rmax" are global min and max values of the radial extent ( input)
! "rint" is the intermediate radius out to which we have uniform mesh.

! "nr_global" is the global number of zones in the radial direction (input).
! "NGhost_global" is the number of ghost cells (input) on the global AMR mesh.
! "nr_geom_zones" is the number of outer zones, which will be log ratioed.

! "rb_save_global" contains boundary information ( "rb_save_global ( 0) = rmin" )
! "dr_save_global" contains width information for each zone.
! "rc_save_global" contains center location for each zone. This is a simple
!  arithmetic averaging of the boundary values. So it is not really a centroid;
!  therefore, do not take it too seriously.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nr_global, NGhost_global, nr_geom_zones, i, ii

      REAL :: rmin, rmax, rint, smallnum, alpha, dr

      REAL
     1 rb_save_global ( 0 - NGhost_global: nr_global + NGhost_global)
      REAL
     1 dr_save_global ( 1  - NGhost_global: nr_global + NGhost_global)
      REAL
     1 rc_save_global ( 1  - NGhost_global: nr_global + NGhost_global)

! ----------------------------------------------------------------------

      dr = ( rint - rmin) / REAL (nr_global - nr_geom_zones)  ! size of uniform inner mesh.

      alpha = LOG ( 1.0 + (rint - rmin)
     1                  / (rint * (nr_global - nr_geom_zones)) )

! This is the suggested value for "rmax"; which is reset in this subroutine.
      rmax = rint * EXP ( alpha * nr_geom_zones)

! Set up "rb_save_global" uniform radial mesh.

      DO i = 0 - NGhost_global, (nr_global - nr_geom_zones)

      rb_save_global ( i) = rmin + dr * REAL ( i)

      END DO

! Set up "rb_save_global" log scaled mesh.

      DO i = (nr_global - nr_geom_zones), nr_global + NGhost_global

      ii = i - (nr_global - nr_geom_zones)

      rb_save_global ( i) = rint * EXP ( alpha * REAL ( ii))

      END DO

! Set up "dr_save_global" and "rc_save_global".

      DO i = 1 - NGhost_global, nr_global + NGhost_global
        dr_save_global ( i) = rb_save_global ( i)
     1                      - rb_save_global ( i - 1)
        rc_save_global ( i) = 0.5 * ( rb_save_global ( i)
     1                                + rb_save_global ( i - 1) )
      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE GENERATE_GLOBAL_RADIAL_UNIF_N_LOG_MESH



