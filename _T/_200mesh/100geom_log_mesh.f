! rm -f x a.F; cp Geom_Log_Mesh_src a.F; gfortran -fdefault-real-8 -o x a.F; ./x



      PROGRAM GEOM_LOG_MESH

      IMPLICIT NONE

! There are "n1" zones in the inner, geometrically ratioed mesh; There are "n2"
! zones in the outer, logarithmetically ratioed mesh.
      INTEGER n1, n2, i, nn

      PARAMETER ( n1 = 100, n2 = 100, nn = 20)

      REAL dr1, ratio_geom, rint, rmax, pow, dr_geom, dr_log

! ----------------------------------------------------------------------

      rmax = 10.0  ! Maximum outer boundary of the mesh.

      DO i = 1, nn

      dr1 = 0.003  ! Size of first zone of the geometrically ratioed mesh.

      ratio_geom = 1.024728 + 0.000005 * (i - 1.0)

! Find outer boundary of the geometrically ratioed mesh and the zone size
! at that boundary.
      rint = 1.0 + dr1 * ( ratio_geom**n1 - 1.0) / ( ratio_geom - 1.0)
      dr_geom = dr1 * ratio_geom**(n1 - 1)

      WRITE ( 6, *)"Using ratio_geom, rint = ", ratio_geom, rint

      pow = 1.0 / REAL ( n2)

      dr_log = rint * ( rmax / rint)**pow - rint

      WRITE ( 6, *)"dr_geom, dr_log = ", dr_geom, dr_log
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------

      STOP

      END PROGRAM GEOM_LOG_MESH

