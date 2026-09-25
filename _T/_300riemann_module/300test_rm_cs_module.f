










      PROGRAM TEST_RM_CS_MODULE

! Written by : Dinshaw Balsara and Sethupathy S.

      USE RIEMANN_3D_CS_AMR_SUPPORT

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, isector, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, ix_patch, iy_patch, iz_patch

      INTEGER nRK_stage, irk, iproblem, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum

      REAL rmin, rmax, time, dtcur, dtcurp1, dtmove, a,
     1     timeimage, dtimage, timestop

!$acc declare create ( time, dtcur, dtcurp1, dtmove )

      REAL sector_vertices ( 4, 3, 6)

! File names used for imaging.

      CHARACTER*8 rhofil, vlxfil, vlyfil, vlzfil, prsfil,
     1            bfxfil, bfyfil, bfzfil, dumpfil,
     1            speciefil

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)::
     1   r4gr


      REAL timing_info1, timing_info2

! Variables used in testing the MPI messaging.

      INTEGER iproc

! Initialize MPI here.

      CALL MPI_INIT (ierr)

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

! Create windows for message passing.
      CALL WINCREATE_MESSAGING ()

! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

! Also create fluid variable arrays that are needed by the algorithm and only need
! to reside on the GPU. These arrays are not directly used for the update
! and will never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,&
!$acc&           xf_area_patch, xf_centroid_patch,                      &

!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &

!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &

!$acc&           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch, &
!$acc&           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch, &
!$acc&           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch, &

!$acc&           v_scale_length,                                        &

!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &

!$acc&           Coef_Basis_Moments,                                    &
!$acc&           itri_patch_LD, itri_patch_RD, itri_patch_RU,           &
!$acc&           itri_patch_LU, istencil_patch,                         &

!$acc&         U, entropygr, Bxf, Byf, Bzf,                             &
!$acc&         V, B1_zc, B2_zc, B3_zc, Wzc, Ex_zc, Ey_zc, Ez_zc,        &

!$acc&         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,       &
!$acc&         divvelgr, msonicgr, flattengr, workgr,                   &
!$acc&         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                   &

!$acc&         V_min_local, V_max_local,                                &

!$acc&         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,          &
!$acc&         EntFlux_x, EntFlux_y, EntFlux_z,                         &
!$acc&         Efx, Efy, Efz)

! Initializing prs and rho floors, inner and outer radius and fluid variables.

      ilevel = 1  ! By default, we work on 1 level.

      OPEN ( UNIT = 6, FILE = "riemann_cs.out")

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

      ix_patch = 2
      iy_patch = 2
      iz_patch = 1

      isector = 1

      rmin = 10.0
      rmax = 15.23

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4     ! Polytropic index for the gas.
      gasconst_euler = 1.0  ! Set up scaled gas constant.
      molewt_euler = 1.4    ! Set up mole wt. of the gas.

      time = 0.0  ! Initial time is always zero.

! Decide which physical problem you want to initialize with "iproblem".
! iproblem = 0 for testing the MPI messaging;
! iproblem = 1 for method of manufactored solution; iproblem = 2 is blast;
! iproblem = 3 is hydro vortex; iproblem = 4 is mhd vortex;
! iproblem = 5 is accuracy test; iproblem = 6 is rotor test;
! iproblem = 7 is Parker wind test; iproblem = 8 hydrodynamic manufactored solution;
! iproblem = 9 is SOD test; iproblem = 10 is Brio-Wu test.

      iproblem = 0


! These parameters control the numerics.

#if ( NUMORDER == 2)
      nRK_stage = 2  ! Two stage RK at second order.
#elif ( NUMORDER == 3)
      nRK_stage = 2  ! Three stage RK at third order. (for now, only 2 stage RK is in code)
#elif ( NUMORDER == 4)
      nRK_stage = 2  ! Four stage RK at fourth order. (for now, only 2 stage RK is in code)
#endif

! Courant number
      cfl_coef = 0.25

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

! Using 100% multid flux can sometimes mean we take flux contributions from
! a stencil that does not depend on the zones that contribute to that face.
! For that reason, we blend multid RS with 1d RS associated with that face.
! Blending in more multidimensionality also seems to help out with pressure
! positivity without any artificial fixes.
! This value should range between [0.4, 0.6] for good performance.

! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
! flux from each edge.
! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
! The electric field is always fully evaluated by the multid. RS.

      blend_multid_flux = 0.1666

! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
! (almost or fully) stationary contacts. To prevent that we blend in a small
! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]

! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
! reverse in smooth flow.

      blend_mdhll_to_mdhllc = 0.0


! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 50000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.
      dtimage = 0.2 / 40.0

! Time at which you want the simulation to stop.
      timestop = 0.2


! By this time we have updated a lot of variables on the CPU. Let us
! migrate that data to the GPU. We still want to continue working
! on CPU, so we set "do_on_gpu = .FALSE.".

      do_on_gpu = .FALSE.

!$acc update device  (/RIEMANN_COM/)


! Set default multispecies values, just to make sure they are not NANS.

      gamma_specie ( :) = gamma
      molewt_specie ( :) = molewt
      heat_of_formation ( :) = 0.0
      cp_specie ( :) = gasconst * gamma_specie ( :)
     1            / ( ( gamma_specie ( :) - 1.0) * molewt_specie ( :) )
      cv_specie ( :) = cp_specie ( :) / gamma_specie ( :)


      gamma_specie_euler ( :) = gamma_euler
      molewt_specie_euler ( :) = molewt_euler
      heat_of_formation_euler ( :) = 0.0
      cp_specie_euler ( :) = gasconst_euler * gamma_specie_euler ( :)
     1 / ( ( gamma_specie_euler ( :) - 1.0) * molewt_specie_euler ( :))
      cv_specie_euler ( :) = cp_specie_euler ( :)
     1                     / gamma_specie_euler ( :)

! Make a default initialization of big fluid arrays. This is still done
! on CPU.

      gamma_eos_gr ( :, :, :) = gamma
      gamma_soundspeed_gr ( :, :, :) = gamma
      mean_mol_wt_gr ( :, :, :) = molewt


! Set the variable names used in the silo file. Not done yet.

! Set file names for the image files. Note that this code just dumps
! out zone-centered variables for imaging. Can't restart from these files.

      rhofil = "rhoa0001"
      vlxfil = "vlxa0001"
      vlyfil = "vlya0001"
      vlzfil = "vlza0001"
      prsfil = "prsa0001"

      bfxfil = "bfxa0001"
      bfyfil = "bfya0001"
      bfzfil = "bfza0001"

      speciefil = "sp010001"

! Set file name for the first dump file.

      dumpfil = "dump0001"


      IF ( me == 1 ) THEN

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Chose physical problem given by iproblem = ",
     1              iproblem
      WRITE ( 6, *)" "


      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"Starting timestep; Documenting constants in the module"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Order of accuracy = ", NUMORDER
      WRITE ( 6, *)"Number of zones theta-phi, Number of ghost cells= ",
     1             ntp_zones, NGhost
      WRITE ( 6, *)"Number of zones in radial direction = ", nr_zones
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Picking one of the patches from the sector"
      WRITE ( 6, *)"isector = ", isector
      WRITE ( 6, *)"ix_patch, iy_patch, iz_patch = ",
     1              ix_patch, iy_patch, iz_patch
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1       "n_cc_components, n_cc_modes = ",
     1        n_cc_components, n_cc_modes
      WRITE ( 6, *)"n_fc_modes = ", n_fc_modes
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"pi, rhofloor, prsfloor, smallnum = ",
     1              pi, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)"gamma, gasconst, molewt = ",gamma, gasconst, molewt
      WRITE ( 6, *)" "
      WRITE ( 6, *)"cfl_coef = ", cfl_coef
      WRITE ( 6, *)" "

      WRITE ( 6, *)"flatten_coef, divratio_einfeldt_rs = ",
     1              flatten_coef, divratio_einfeldt_rs
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lapidus_threshold, lapidus_coef = ",
     1              lapidus_threshold, lapidus_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"blend_multid_flux, blend_mdhll_to_mdhllc = ",
     1              blend_multid_flux, blend_mdhll_to_mdhllc

      WRITE ( 6, *)"ntstep, ndumpstep, nimgstep, time, timestop = ",
     1              ntstep, ndumpstep, nimgstep, time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of RK stages, nRK_stage = ", nRK_stage
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

      END IF ! End of "IF ( me == 1 ) THEN".


! Here we build sibling info. for the base level mesh.

      CALL BUILD_BASE_PATCHES ()

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! Write out the iactive, isibling_* data.

!      CALL TEST_1_PRINT_IACTIVE_ISIBLING ( )

! Build a 1D mesh in the radial direction. This calls a generic on-processor
! serial code. Therefore, synchronize after exit.

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes,

     1           ii_mode, jj_mode, kk_mode)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)






! Loop over patches, building their geometry and initializing MHD
! variables for all patches on this processor.

! Since we are dealing with only one level mesh, each processor will work on
! "npatch_per_proc" number of patches. This variable has already been built.
! We loop over patches here.

! Note that "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! and "iactive_FromProc" to be properly initialized on each processor. This
! was already done via calll to "BUILD_BASE_PATCHES".


      DO ii = 1, npatch_per_proc


      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

! From "ipatch" (which is the logical patch number in "isector") we can
! extract "ix_patch, iy_patch, iz_patch".

      jj = ipatch
      iz_patch = 1 + (jj - 1) / (ntp_patches**2)

      jj = jj - (iz_patch - 1) * ntp_patches**2
      iy_patch = 1 + (jj - 1) / ntp_patches

      jj = jj - (iy_patch - 1) * ntp_patches
      ix_patch = jj

! Fill in "isector" and level information. Also fill in "ix_patch, iy_patch, iz_patch"
! information.
      isector_ilevel_amr ( 1, ii) = isector
      isector_ilevel_amr ( 2, ii) = 1  ! For now, only one level in the AMR hierarchy.

! Fill in "ix_patch, iy_patch, iz_patch" for this patch
      ixyz_patch_amr ( 1, ii) = ix_patch
      ixyz_patch_amr ( 2, ii) = iy_patch
      ixyz_patch_amr ( 3, ii) = iz_patch

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! Save "vert_coords_patch, imask_active_patch" for later reuse. These are
! also needed for imaging.

      vert_coords_amr ( :, :, :, :, ii) =
     1                              vert_coords_patch ( :, :, :, :)

      imask_active_amr ( :, :, :, ii) = imask_active_patch ( :, :, :)


      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Understand that we are working on a single patch on each processor.
! We have just filled in "vert_coords_patch, imask_active_patch" which is
! the minimal amount of information needed for regenating geometry,
! quadrature points, weights etc.

! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! Save information on whether corners are tri-corners or quad-corners.
! Also, save stencil information because for cube-sphere this is only 2D
! and very lightweight.

      itri_patch_amr ( 1, ii) = itri_patch_LD
      itri_patch_amr ( 2, ii) = itri_patch_RD
      itri_patch_amr ( 3, ii) = itri_patch_RU
      itri_patch_amr ( 4, ii) = itri_patch_LU

      istencil_amr ( :, :, :, :, ii) = istencil_patch ( :, :, :, :)


      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".

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


! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)


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

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! By this point, all geometry is built. Therefore, we are now in a position
! to initialize the MHD variables for this patch.

      a = cfl_coef
      cfl_coef = 0.05

!$acc update device ( cfl_coef)

      CALL TEST_MPI_INIT_PATCH_CUBE_SPHERE

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

     1           V_min_local, V_max_local)

      cfl_coef = a  ! Reset the CFL number to its original value.

! By this point, all "U, V, entropygr, Bxf, Byf, Bzf, gamma_eos_gr" for this
! "ipatch" are initialized in the interior of the mesh. But each processor
! can have multiple patches, which are all stored in
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

      Fluct_x_Pl_amr ( :, :, :, :, ii) = Fluct_x_Pl ( :, :, :, :)
      Fluct_y_Pl_amr ( :, :, :, :, ii) = Fluct_y_Pl ( :, :, :, :)
      Fluct_z_Pl_amr ( :, :, :, :, ii) = Fluct_z_Pl ( :, :, :, :)

      Fluct_x_Mi_amr ( :, :, :, :, ii) = Fluct_x_Mi ( :, :, :, :)
      Fluct_y_Mi_amr ( :, :, :, :, ii) = Fluct_y_Mi ( :, :, :, :)
      Fluct_z_Mi_amr ( :, :, :, :, ii) = Fluct_z_Mi ( :, :, :, :)

      EntFlux_x_amr ( :, :, :, ii) = EntFlux_x ( :, :, :)
      EntFlux_y_amr ( :, :, :, ii) = EntFlux_y ( :, :, :)
      EntFlux_z_amr ( :, :, :, ii) = EntFlux_z ( :, :, :)

      Efx_amr ( :, :, :, ii) = Efx ( :, :, :)
      Efy_amr ( :, :, :, ii) = Efy ( :, :, :)
      Efz_amr ( :, :, :, ii) = Efz ( :, :, :)


      END DO  ! End of "DO ii = 1, npatch_per_proc".

! End of looping over patches, building their geometry and initializing MHD
! variables for all patches on this processor.

! Test the XY messaging
!       CALL TEST_XY_MESSAGING ()

! Calculate error in Ghost zones before MPI messaging.

      iproc = 1

      iproc = MAX ( iproc, 1)
      iproc = MIN ( iproc, nimages)

      IF ( me == iproc ) THEN
        WRITE ( 6, *)"iproc, nimages, npatch_per_proc = ",
     1                iproc, nimages, npatch_per_proc
        WRITE ( 6, *)"----------------------------------------------"
        WRITE ( 6, *)"Before MPI messaging:"
        WRITE ( 6, *)" "
      END IF

      CALL TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc,
     1                 time, dtcur, rmin, rmax )

! Do messaging and fill in Ghost zones.

! Fill ghoszt zones and faces in X and Y directions; i.e., along the
! spherical surface.

! Fill send buffers in both X and Y directions.
      CALL LOAD_LOCAL_BUFFERS_UXY ()

! Fetch the X-data from left and right patches.
      CALL FETCH_BUFFERS_FILL_IN_UX ()

! Fill send buffers in both X and Y directions.
      CALL LOAD_LOCAL_BUFFERS_UXY ()

! Fetch the Y-data from up and down patches.
      CALL FETCH_BUFFERS_FILL_IN_UY ()

! Till above all the ghost zones are filled for interior and the bottom
! patches in the sector. For the top patches, we need to call the filling,
! fecthing subroutines in Y-direction again.

      CALL LOAD_LOCAL_BUFFERS_UXY ()

      CALL FETCH_BUFFERS_FILL_IN_UY_TOP ()

! Now message in the Z direction. This completes the messaging.

! Fill send buffers in Z direction.
      CALL LOAD_LOCAL_BUFFERS_UZ ()

! Fetch the Z-data from top and bottom patches.
      CALL FETCH_BUFFERS_FILL_IN_UZ ()

! Synchronize Fluxes and E-fields in the same sequence as of MHD variables
! written above.

      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()

      CALL FETCH_BUFFERS_FILL_IN_Flux_X ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Y ()

      CALL LOAD_LOCAL_BUFFERS_Flux_Z ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Z ()

      CALL SYNCHRONIZE_EFZ_QUAD_CORNER ()

! Assert the ghost zones and faces have received the correct values
! for the case of "iproblem = 0".

      IF ( me == iproc ) THEN
        WRITE ( 6, *)"----------------------------------------------"
        WRITE ( 6, *)"After MPI messaging: "
        WRITE ( 6, *)" "
      END IF

      CALL TEST_MESSAGING_ERROR_CALCULATION ( iproblem, iproc,
     1                 time, dtcur, rmin, rmax )

! All the ghost zones and faces are filled. Now free up all the handles.
      CALL WINFREE_MESSAGING ()

! Set the variable names used in the silo file.

      silo_allvars ( 1) = "Density "
      silo_allvars ( 2) = "XMomentum "
      silo_allvars ( 3) = "YMomentum "
      silo_allvars ( 4) = "ZMomentum "
      silo_allvars ( 5) = "Energy "
      silo_allvars ( 6) = "XMagneticField "
      silo_allvars ( 7) = "YMagneticField "
      silo_allvars ( 8) = "ZMagneticField "

#if( IPDE_SYSTEM == 1 && NFLUID >= 1)
      silo_allvars ( 8 + 1) = "Tracer_1 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 2)
      silo_allvars ( 8 + 2) = "Tracer_2 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 3)
      silo_allvars ( 8 + 3) = "Tracer_3 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 4)
      silo_allvars ( 8 + 4) = "Tracer_4 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 5)
      silo_allvars ( 8 + 5) = "Tracer_5 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 6)
      silo_allvars ( 8 + 6) = "Tracer_6 "
#endif

! ----------

#if( IPDE_SYSTEM == 2)

      silo_allvars ( 9) = "Pdiff_Ppar_Pperp "

#if ( NVAR >= 10)
      silo_allvars ( 10) = "Tracer_1 "
#endif

#if ( NVAR >= 11)
      silo_allvars ( 11) = "Tracer_2 "
#endif

#endif
! End of "IPDE_SYSTEM == 2".

! ----------

#if( IPDE_SYSTEM == 3)

      silo_allvars ( 9) = "Pdiff_Ppar_Pperp "

      silo_allvars ( 10) = "Prs_electron "

#if ( NVAR >= 11)
      silo_allvars ( 11) = "Tracer_1 "
#endif

#if ( NVAR >= 12)
      silo_allvars ( 12) = "Tracer_2 "
#endif

#if ( NVAR >= 13)
      silo_allvars ( 13) = "Tracer_3 "
#endif

#if ( NVAR >= 14)
      silo_allvars ( 14) = "Tracer_4 "
#endif

#if ( NVAR >= 15)
      silo_allvars ( 15) = "Tracer_5 "
#endif

#if ( NVAR >= 16)
      silo_allvars ( 16) = "Tracer_6 "
#endif


#endif
! End of "IPDE_SYSTEM == 3".

! ----------
! ----------

#if( USE_B0B1_FORMULATION == 1)
      silo_allvars ( NVAR + NFLUID + 1) = "XTotalMagnField "
      silo_allvars ( NVAR + NFLUID + 2) = "YTotalMagnField "
      silo_allvars ( NVAR + NFLUID + 3) = "ZTotalMagnField "
#endif

! ----------------------------------------
! Visualizing the extra internal variables.

! --------------------
#if( USE_B0B1_FORMULATION == 1)
!-----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NVAR + NFLUID + 3 + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NVAR + NFLUID + 3 + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NVAR + NFLUID + 3 + 3) = "xflag_entflag "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NVAR + NFLUID + 3 + 4) = "xflag_tempfloor "
#endif

! ----------
#else
! ----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NVAR + NFLUID + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NVAR + NFLUID + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NVAR + NFLUID + 3) = "xflag_entprs "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NVAR + NFLUID + 4) = "xflag_tempfloor "
#endif

! ----------
#endif
! --------------------


#if ( SILO_EXTRA_DATA >= 5)
        WRITE( 6, *)" "
        WRITE( 6, *)"Stopping at - Set the silo variable names"
        WRITE( 6, *)"SILO_EXTRA_DATA is expected to be <= 4"
        WRITE( 6, *)"SILO_EXTRA_DATA = ", SILO_EXTRA_DATA
        WRITE( 6, *)" "
        CALL FLUSH( 6)
        ERROR STOP
#endif

! ----------------------------------------------------------------------

      is_bottom_or_top ( :, :) = 0.0

      DO ii = 1, npatch_per_proc

        iz_patch = ixyz_patch_amr ( 3, ii)

        IF ( iz_patch == 1)           is_bottom_or_top ( 1, ii) = 1
        IF ( iz_patch == nr_patches)  is_bottom_or_top ( 2, ii) = 1

      END DO

      filenum = 1;  time = 0.0

      CALL  WRITEDATA ( npatch_per_proc,
     1                  iactive_FromProc ( :, 2),
     1                  nr_zones, ntp_zones, NGhost,
     1                  is_bottom_or_top,
     1                  itri_patch_amr,
     1                  vert_coords_amr,
     1                  n_cc_components,
     1                  SILO_VARLEN, silo_allvars,
     1                  U_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                  filenum, time
     1                )


      CALL MPI_FINALIZE (ierr)


      END PROGRAM TEST_RM_CS_MODULE

