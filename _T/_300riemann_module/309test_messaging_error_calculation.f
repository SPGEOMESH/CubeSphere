








      SUBROUTINE TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc,
     1                 time, dtcur, rmin, rmax )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.
!
! This test subroutine computes the error in "U_amr, entropygr_amr,
! Bxf_amr, Byf_amr, Bzf_amr".
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iproblem, iproc

      REAL time, dtcur, rmin, rmax

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ix_patch, iy_patch, iz_patch, isector,
     1        itri_patch_LD, itri_patch_LU,
     1        itri_patch_RD, itri_patch_RU

      INTEGER ixx, iyy, izz, i, ii, izb, ize, im1, ip1

      REAL xc, yc, zc

      REAL err_U ( n_cc_components),
     1     err_Bxf, err_Byf, err_Bzf,
     1     err_Fx, err_Fy, err_Fz,
     1     err_EntFx, err_EntFy, err_EntFz,
     1     err_Efx, err_Efy, err_Efz,

     1     err_Flx_Pl, err_Fly_Pl, err_Flz_Pl,
     1     err_Flx_Mi, err_Fly_Mi, err_Flz_Mi,

! Global error from all processors.
     1     gl_err_U ( n_cc_components),
     1     gl_err_Bxf, gl_err_Byf, gl_err_Bzf,
     1     gl_err_Fx, gl_err_Fy, gl_err_Fz,
     1     gl_err_EntFx, gl_err_EntFy, gl_err_EntFz,
     1     gl_err_Efx, gl_err_Efy, gl_err_Efz,

     1     gl_err_Flx_Pl, gl_err_Fly_Pl, gl_err_Flz_Pl,
     1     gl_err_Flx_Mi, gl_err_Fly_Mi, gl_err_Flz_Mi

      INTEGER :: write_data = 0

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------

      IF ( me == 1 ) THEN

      IF ( iproblem .NE. 0 ) THEN
        WRITE ( 6, *)"Stopping at TEST_MESSAGING_ERROR_CALCULATION"
        WRITE ( 6, *)"Expecting iproblem == 0"
        WRITE ( 6, *)"iproblem = ",iproblem
        STOP
      END IF

      END IF

! Initialize the error variables to zero.
      err_U ( :) = 0.0;   gl_err_U ( :) = 0.0

      err_Bxf = 0.0;      err_Byf = 0.0;      err_Bzf = 0.0
      gl_err_Bxf = 0.0;   gl_err_Byf = 0.0;   gl_err_Bzf = 0.0

      err_Fx = 0.0;      err_Fy = 0.0;     err_Fz = 0.0
      err_EntFx = 0.0;   err_EntFy = 0.0;  err_EntFz = 0.0
      err_Efx = 0.0;     err_Efy = 0.0;    err_Efz = 0.0

      gl_err_Fx = 0.0;      gl_err_Fy = 0.0;     gl_err_Fz = 0.0
      gl_err_EntFx = 0.0;   gl_err_EntFy = 0.0;  gl_err_EntFz = 0.0
      gl_err_Efx = 0.0;     gl_err_Efy = 0.0;    gl_err_Efz = 0.0

      err_Flx_Pl = 0.0;     err_Fly_Pl = 0.0;     err_Flz_Pl = 0.0
      err_Flx_Mi = 0.0;     err_Fly_Mi = 0.0;     err_Flz_Mi = 0.0

      gl_err_Flx_Pl = 0.0;  gl_err_Fly_Pl = 0.0;  gl_err_Flz_Pl = 0.0
      gl_err_Flx_Mi = 0.0;  gl_err_Fly_Mi = 0.0;  gl_err_Flz_Mi = 0.0

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! sector number and patch location in a sector.

      isector  = isector_ilevel_amr ( 1, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! ----------------------------------------------------------------------

      vert_coords_patch( :, :, :, :) = vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) = imask_active_amr ( :, :, :, ii)

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".
! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Build the quadrature points, weights and normals for the integration
! of magnetic field in the Multi-D Riemann solver.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Compute error in "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------
! Change the ghost zone inclusion along "z-axis", if the patch is at inner or
! outer r-boundary of the simulation, where there is no MPI messaging.
      izb = 1 - NGhost
      ize = nr_zones + NGhost

      IF ( iz_patch == 1         ) izb = 1         ! patch at inner r-boundary.
      IF ( iz_patch == nr_patches) ize = nr_zones  ! patch at outer r-boundary.

! ----------------------------------------------------------------------

!      IF ( me == iproc ) write_data = 1

      U ( :, :, :, :, 1) = U_amr ( :, :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)

      Fluct_x_Pl ( :, :, :, :) = Fluct_x_Pl_amr ( :, :, :, :, ii)
      Fluct_y_Pl ( :, :, :, :) = Fluct_y_Pl_amr ( :, :, :, :, ii)
      Fluct_z_Pl ( :, :, :, :) = Fluct_z_Pl_amr ( :, :, :, :, ii)

      Fluct_x_Mi ( :, :, :, :) = Fluct_x_Mi_amr ( :, :, :, :, ii)
      Fluct_y_Mi ( :, :, :, :) = Fluct_y_Mi_amr ( :, :, :, :, ii)
      Fluct_z_Mi ( :, :, :, :) = Fluct_z_Mi_amr ( :, :, :, :, ii)

      EntFlux_x ( :, :, :) = EntFlux_x_amr ( :, :, :, ii)
      EntFlux_y ( :, :, :) = EntFlux_y_amr ( :, :, :, ii)
      EntFlux_z ( :, :, :) = EntFlux_z_amr ( :, :, :, ii)

      Efx ( :, :, :) = Efx_amr ( :, :, :, ii)
      Efy ( :, :, :) = Efy_amr ( :, :, :, ii)
      Efz ( :, :, :) = Efz_amr ( :, :, :, ii)

      CALL TEST_MPI_ERROR_IN_B_XYZ_FACES

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_quadpts_patch, v_quadwt_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           EntFlux_x, EntFlux_y, EntFlux_z, Efx, Efy, Efz,

     1           Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1           Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi,

     1           gamma_eos_gr,

     1           izb, ize, write_data,

     1           err_U, err_Bxf, err_Byf, err_Bzf,
     1           err_Fx, err_Fy, err_Fz,
     1           err_EntFx, err_EntFy, err_EntFz,
     1           err_Efx, err_Efy, err_Efz,

     1           err_Flx_Pl, err_Fly_Pl, err_Flz_Pl,
     1           err_Flx_Mi, err_Fly_Mi, err_Flz_Mi)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_U" from all processors.
! Then broadcast the global maximum "gl_err_U" to all processors from
! my_rank == 0  processor.

      CALL MPI_REDUCE ( err_U, gl_err_U, n_cc_components, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast   ( gl_err_U, n_cc_components, MPI_REAL8, 0,
     1                   MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_Bxf, err_Byf, err_Bzf" from all processors.
! Then broadcast the global maximum "gl_err_Bxf, gl_err_Byf, gl_err_Bzf"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Bxf, gl_err_Bxf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Byf, gl_err_Byf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Bzf, gl_err_Bzf, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_Bxf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Byf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Bzf, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_Fx, err_Fy, err_Fz" from all processors.
! Then broadcast the global maximum "gl_err_Fx, gl_err_Fy, gl_err_Fz"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Fx, gl_err_Fx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fy, gl_err_Fy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fz, gl_err_Fz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_Fx, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fy, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fz, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_Flx_Pl, err_Fly_Pl, err_Flz_Pl" from all procs.
! Then broadcast the global maximum "gl_err_Flx_Pl, gl_err_Fly_Pl, gl_err_Flz_Pl"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Flx_Pl, gl_err_Flx_Pl, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fly_Pl, gl_err_Fly_Pl, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Flz_Pl, gl_err_Flz_Pl, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_Flx_Pl, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fly_Pl, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Flz_Pl, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

! ----------------------------------------
! Reduce and get max error for "err_Flx_Mi, err_Fly_Mi, err_Flz_Mi" from all procs.
! Then broadcast the global maximum "gl_err_Flx_Mi, gl_err_Fly_Mi, gl_err_Flz_Mi"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Flx_Mi, gl_err_Flx_Mi, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Fly_Mi, gl_err_Fly_Mi, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Flz_Mi, gl_err_Flz_Mi, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_Flx_Mi, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Fly_Mi, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Flz_Mi, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_EntFx, err_EntFy, err_EntFz" from all processors.
! Then broadcast the global maximum "gl_err_EntFx, gl_err_EntFy, gl_err_EntFz"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_EntFx, gl_err_EntFx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_EntFy, gl_err_EntFy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_EntFz, gl_err_EntFz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_EntFx, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_EntFy, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_EntFz, 1, MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Reduce and get max error for "err_Efx, err_Efy, err_Efz" from all processors.
! Then broadcast the global maximum "gl_err_Efx, gl_err_Efy, gl_err_Efz"
! to all processors from my_rank == 0  processor.

      CALL MPI_REDUCE ( err_Efx, gl_err_Efx, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Efy, gl_err_Efy, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE ( err_Efz, gl_err_Efz, 1, MPI_REAL8,
     1                  MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast ( gl_err_Efx, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Efy, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast ( gl_err_Efz, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

      IF ( me == iproc ) THEN

      WRITE ( 6, *)"err_U ( 2) = ",err_U ( 2)
      WRITE ( 6, *)"err_U ( 3) = ",err_U ( 3)
      WRITE ( 6, *)"err_U ( 4) = ",err_U ( 4)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_U ( 2) = ",gl_err_U ( 2)
      WRITE ( 6, *)"gl_err_U ( 3) = ",gl_err_U ( 3)
      WRITE ( 6, *)"gl_err_U ( 4) = ",gl_err_U ( 4)
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Bxf = ",err_Bxf
      WRITE ( 6, *)"err_Byf = ",err_Byf
      WRITE ( 6, *)"err_Bzf = ",err_Bzf
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Bxf = ",gl_err_Bxf
      WRITE ( 6, *)"gl_err_Byf = ",gl_err_Byf
      WRITE ( 6, *)"gl_err_Bzf = ",gl_err_Bzf
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Fx = ",err_Fx
      WRITE ( 6, *)"err_Fy = ",err_Fy
      WRITE ( 6, *)"err_Fz = ",err_Fz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Fx = ",gl_err_Fx
      WRITE ( 6, *)"gl_err_Fy = ",gl_err_Fy
      WRITE ( 6, *)"gl_err_Fz = ",gl_err_Fz
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Flx_Pl = ",err_Flx_Pl
      WRITE ( 6, *)"err_Fly_Pl = ",err_Fly_Pl
      WRITE ( 6, *)"err_Flz_Pl = ",err_Flz_Pl
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Flx_Pl = ",gl_err_Flx_Pl
      WRITE ( 6, *)"gl_err_Fly_Pl = ",gl_err_Fly_Pl
      WRITE ( 6, *)"gl_err_Flz_Pl = ",gl_err_Flz_Pl
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Flx_Mi = ",err_Flx_Mi
      WRITE ( 6, *)"err_Fly_Mi = ",err_Fly_Mi
      WRITE ( 6, *)"err_Flz_Mi = ",err_Flz_Mi
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Flx_Mi = ",gl_err_Flx_Mi
      WRITE ( 6, *)"gl_err_Fly_Mi = ",gl_err_Fly_Mi
      WRITE ( 6, *)"gl_err_Flz_Mi = ",gl_err_Flz_Mi
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_EntFx = ",err_EntFx
      WRITE ( 6, *)"err_EntFy = ",err_EntFy
      WRITE ( 6, *)"err_EntFz = ",err_EntFz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_EntFx = ",gl_err_EntFx
      WRITE ( 6, *)"gl_err_EntFy = ",gl_err_EntFy
      WRITE ( 6, *)"gl_err_EntFz = ",gl_err_EntFz
      WRITE ( 6, *)"------------------------------------------------"
      WRITE ( 6, *)"err_Efx = ",err_Efx
      WRITE ( 6, *)"err_Efy = ",err_Efy
      WRITE ( 6, *)"err_Efz = ",err_Efz
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gl_err_Efx = ",gl_err_Efx
      WRITE ( 6, *)"gl_err_Efy = ",gl_err_Efy
      WRITE ( 6, *)"gl_err_Efz = ",gl_err_Efz
      WRITE ( 6, *)"------------------------------------------------"

      END IF

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

      END SUBROUTINE TEST_MESSAGING_ERROR_CALCULATION

