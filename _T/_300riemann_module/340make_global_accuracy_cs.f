







      SUBROUTINE MAKE_GLOBAL_ACCURACY_CS ( iproblem,
     1                                     rmin, rmax, time, dtcur )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.
!
! Compute the accuracy of the simulation. It prints out the error after
! each time step.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER iproblem

      REAL rmin, rmax, time, dtcur

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, isector, ix_patch, iy_patch, iz_patch

      INTEGER nRK_stage, irk, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum, ibdy_inner, ibdy_outer

      INTEGER, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost, 2) :: bcarr_r

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Now we loop over all "npatch_per_proc" patches on each processor and
! build fluxes and electric fields.

      V_min_reduce ( :) = 1.0e32
      V_max_reduce ( :) = -1.0e32

      U_L1_reduce      ( :) = 0.0
      U_Linf_reduce    ( :) = 0.0
      U_zc_Linf_reduce ( :) = 0.0

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! STEP 1) First reverse all that we wrote into the big "*_amr" arrays during
! initialization.

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! ------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! ------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! ------------------------------

      U ( :, :, :, 1:n_cc_components, 1) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr ( :, :, :, 1) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)

      flattengr ( :, :, :) = flatten_amr ( :, :, :, ii)

! ----------------------------------------------------------------------

!$acc update device (

!$acc&       time, dtcur, cfl_coef,                                     &

!$acc&       isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch,     &
!$acc&       vert_coords_patch, imask_active_patch,                     &

!$acc&       itri_patch_LD, itri_patch_RD, itri_patch_RU,               &
!$acc&       itri_patch_LU, istencil_patch,                             &

!$acc&       U, entropygr, Bxf, Byf, Bzf, flattengr )


! ----------------------------------------------------------------------
! STEP 2) Now build all the rest of the geometry information for this patch.

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the area-related quadrature points,
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

! Build the quadrature points, weights and normals for the integration
! of magnetic field in the Multi-D Riemann solver.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! --------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! --------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! STEP 3) With geometry built, now do limiting, predictor step, flux
! and field evaluation.

! --------------------------------------------------

! Pad the boundaries in the radial direction if this is the innermost
! radial patch or the outermost radial patch.

! ------------------------------

      IF ( iz_patch .EQ. 1) THEN

! --------------------
! Inner radial boundary needs to the padded.

! Most of the geomesh simulations require the re-initialization of boundary
! conditions at the inner boundary. The ADD_BOUNDARY subroutine does that.

      ibdy_inner = 1   ! 1 == Impose Add_boundary at inner; 0 == No Add_Bdy
      ibdy_outer = 1   ! 1 == Impose Add_boundary at outer; 0 == No Add_Bdy

!$acc update device ( ibdy_inner, ibdy_outer)

      CALL ADD_BOUNDARY

     1         ( ibdy_inner, ibdy_outer,

     1           nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_centroid_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           gamma_eos_gr, bcarr_r,

     1           V_min_local, V_max_local)

! --------------------

      END IF

! ------------------------------

      IF ( iz_patch .EQ. nr_patches) THEN
! Outer radial boundary needs to the padded.

! Most of the geomesh simulations require to implement the outflow  boundary
! conditions at the outer boundary. The PAD_BOUNDARY subroutine does that.

      bcarr_r ( :, :, 1)  = 17   ! No inner boundary enforcement here.
      bcarr_r ( :, :, 2)  = 17    ! Enforce outflow conditions at Outer boundary.

!$acc update device ( bcarr_r)

      CALL PAD_BOUNDARY_R

     1    ( nr_zones, ntp_zones, NGhost,
     1      n_cc_components, n_cc_modes,
     1      n_fc_modes,

     1      vert_coords_patch, imask_active_patch,

     1      xf_area_patch, yf_area_patch, zf_area_patch,
     1      v_volume_patch,

     1      v_centroid_patch,

     1      U, entropygr, Bxf, Byf, Bzf,

     1      bcarr_r)

      END IF

! --------------------------------------------------

! Apply limiter to MHD variables.

#if (NUMORDER == 2)
      CALL LIMIT_ALL_SLOPES_P1
#elif (NUMORDER == 3)
      CALL LIMIT_ALL_SLOPES_P12
#endif
     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iz_patch, nr_patches,

     1         ii_mode, jj_mode, kk_mode,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc, Wzc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! Store local min and max for each patch that resides on this processor.
      V_min_amr ( :, ii) = V_min_local ( :)
      V_max_amr ( :, ii) = V_max_local ( :)

! For all the patches that reside on this processor, make a reduction of their
! min and max values.
      V_min_reduce = AMIN1 ( V_min_reduce, V_min_local)
      V_max_reduce = AMAX1 ( V_max_reduce, V_max_local)

! --------------------------------------------------
! Compute the error.

      CALL ACCURACY_CALCULATION_PATCH_CS

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           ii_mode, jj_mode, kk_mode,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_quadpts_patch, v_quadwt_patch,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           EntFlux_x, EntFlux_y, EntFlux_z, Efx, Efy, Efz,
     1           gamma_eos_gr,

     1           U_L1, U_Linf, U_zc_Linf)

      U_L1_reduce ( :) = U_L1_reduce ( :) + U_L1 ( :)

      U_Linf_reduce ( :) = MAX ( U_Linf_reduce ( :), U_Linf ( :))

      U_zc_Linf_reduce ( :) = MAX ( U_zc_Linf_reduce( :), U_zc_Linf( :))

! --------------------
#if 0
      IF ( me == 1) THEN

        WRITE ( 6, *)"U_L1       = ", U_L1
        WRITE ( 6, *)"U_Linf     = ", U_Linf
        WRITE ( 6, *)"U_zc_Linf  = ", U_zc_Linf

      END IF
#endif
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      U_L1_reduce ( :) = U_L1_reduce ( :) / REAL ( npatch_per_proc)

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------

! Now make a reduction across all processors for "V_min_reduce, V_max_reduce".

      U_L1      ( :) = U_L1_reduce      ( :)
      U_Linf    ( :) = U_Linf_reduce    ( :)
      U_zc_Linf ( :) = U_zc_Linf_reduce ( :)


      CALL MPI_ALLREDUCE ( U_L1,              ! Input for reudction
     1                     U_L1_reduce,       ! Global sum
     1                     n_cc_components,   ! length of reduction variable
     1                     MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      U_L1_reduce ( :) = U_L1_reduce ( :) / REAL ( nimages)


      CALL MPI_ALLREDUCE ( U_Linf,            ! Input for reudction
     1                     U_Linf_reduce,     ! Global maximum
     1                     n_cc_components,   ! length of reduction variable
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE ( U_zc_Linf,         ! Input for reudction
     1                     U_zc_Linf_reduce,  ! Global maximum
     1                     n_cc_components,   ! length of reduction variable
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)


! Write out the global min and max of all variables from processor 0.

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" ---------------------- Error --------------- "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_L1      = ", U_L1_reduce ( :)
      WRITE ( 6, *)"U_Linf    = ", U_Linf_reduce ( :)
      WRITE ( 6, *)"U_zc_Linf = ", U_zc_Linf_reduce ( :)
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE MAKE_GLOBAL_ACCURACY_CS


