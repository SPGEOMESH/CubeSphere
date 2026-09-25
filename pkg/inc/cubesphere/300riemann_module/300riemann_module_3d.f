








      MODULE RIEMANN_3D_CS_AMR_SUPPORT

! ----------------------------------------------------------------------
! Written by : Dinshaw Balsara

! 3D MPI-3 and Fortran-based Adaptive Mesh Refinement version of Cubed Sphere
! based MHD code from Dinshaw Balsara's group. Uses Silo for I/O.
! An MPI-twinned Fortran compiler, along with Silo libraries, is all one needs to run.

! Sectors and Patches:-
! Code operates on six spherical sectors, each of which may be split into
! many (logically) cubical patches. Each patch has "ntp_zonesXntp_zones"
! in the theta-phi directions and "nr_zones" in the radial direction.
! Each patch has "NGhost" number of ghost zones.

! Each Sector:-
! Within each sector, the base-level mesh can have "ntp_patches X ntp_patches" in
! the theta-phi directions and "nr_patches" in the radial direction. "ilevel_max"
! number of levels can be present; so that each base level zone can be split
! into up to "iscale**3" number of sub-cells.

! Global Mesh:-
! The above information can be used to build "ntp_global", the number of zones
! in any one theta-phi direction on any of the 6 sectors. Similarly, we have
! "nr_global" and "NGhost_global". Each sector has a complete 2D lattice
! on the unit sphere at the finest level, so that the spherical geometry can
! be extracted at any desired level from that.

! On the base mesh, we also retain "ntp_global_base" and "nr_global_base" so
! that building sibling relationships is very easy for one-level mesh.

! Active Patches:-
! "Max_base_patches_per_sector" is the max patches per sector on the base level.

! "Max_Theoretical_patches_per_sector" is the max number of theoretical patches
! that we can have in one sector on a "ilevel_max" global hierarchy.
! We keep track of whether that patch location is filled or not. If it
! is filled, we provide its patch number in that sector.

! Reality is that we only have "Max_patches_per_proc" on each processor.
! This is V. V. IMP.:- "Max_patches_per_proc" is the number of patches per core.
! Since we don't know how many cores the user will be running on, you have to set it
! by hand before each run.

! Processor numbering starts from 1; so remember to subtract 1 for MPI's C format.

! MHD Variables:-
! Thre are "n_cc_components" number of zone-centered variables which are
! reconstructed (based on order) to have up to "n_cc_modes". Each face has number of
! vector field "n_fc_components" which can each be reconstructed to have
! "n_fc_modes". After div-preserving reconstruction, we will have "n_cc_modes"
! number of modes for each component of the vector field.

! Parallelism:-
! The maximum number of cores is just "6 * ntp_patches**2 * nr_patches". We can also put
! 2 patches on each core with "6 * ntp_patches**2 * nr_patches / 2" cores or 3 patches
! per core with "6 * ntp_patches**2 * nr_patches / 3" cores and so on.
! Only thing is: Once you have chosen a number of cores, you can't change it.

! Overall zones in each sector:-
! The overall zones in each sector will be "ntp_patches * ntp_zones" in each of
! the angular, i.e. theta-phi, directions. In the radial direction we have
! "nr_patches * nr_zones" for a total of
! "(ntp_patches * ntp_zones)**2 * (nr_patches * nr_zones)" zones per sector.

! ----------------------------------------------------------------------


      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT NONE

      include 'mpif.h'

! List of public routines for this module:

      PUBLIC :: BUILD_BASE_PATCHES,

! Subroutines for creating MPI Windows and freeing them up.
     1      WINCREATE_MESSAGING, WINFREE_MESSAGING,

! Load Local Buffers; first load them in x and y-directions. Then load in z-direction.
     1      LOAD_LOCAL_BUFFERS_UXY, LOAD_LOCAL_BUFFERS_UZ,

     1      LOAD_LOCAL_BUFFERS_Flux_XY, LOAD_LOCAL_BUFFERS_Flux_Z,

! Fetch buffers from neighbors and initialize suitably.
     1      FETCH_BUFFERS_FILL_IN_UX, FETCH_BUFFERS_FILL_IN_UY,
     1      FETCH_BUFFERS_FILL_IN_UZ, FETCH_BUFFERS_FILL_IN_UY_TOP,

     1      FETCH_BUFFERS_FILL_IN_Flux_X, FETCH_BUFFERS_FILL_IN_Flux_Y,
     1      FETCH_BUFFERS_FILL_IN_Flux_Z,
     1      SYNCHRONIZE_EFZ_QUAD_CORNER,

     1      LOAD_AMR_FLUXES_AT_BDY

! todo NB I moved these out of the module, so it might be required to munge the tests before they work.
! Testing Subroutine, Print amr variable values in a sector.
!    1      PRINT_VARIABLES_ONE_PATCH_PER_SECTOR,
!    1      TEST_XY_MESSAGING, TEST_MESSAGING_ERROR_CALCULATION,
!!if( ACCURACY_ANALYSIS == 1 )
!    1      MAKE_GLOBAL_ACCURACY_CS,
!endif
!! Silo subroutines.
!    1      WRITE_SILO_ONE_BLOCK

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! START LIST of parameters that control Cube Sphere AMR mesh and its processing.

! ---------------
! Variables
! ---------------


      INTEGER, BIND(C) :: ntp_zones, nr_zones, NGhost, ntp_patches, nr_patches,

     1        ilevel_max, iscale,

     1        ntp_global, nr_global, NGhost_global,

     1        ntp_global_base, nr_global_base,

     1        Max_base_patches_per_sector,

     1        Max_Theoretical_patches_per_sector, Max_patches_per_proc,

     1        n_cc_components, n_cc_modes,
     1        n_fc_components, n_fc_modes



      PARAMETER

! Variables that specify the number of zones "ntp_zones, nr_zones" and the
! number of ghost zones "NGhost" for a single patch.
     1          ( ntp_zones = 35, nr_zones = 45, NGhost = 4,

! Number of patches "ntp_patchesXntp_patchesXnr_patches" within a sector on
! the global base-level mesh.
     1            ntp_patches = 2, nr_patches = 16,

! Maximum level in AMR hierarchy and scale between base and finest meshes.
     1            ilevel_max = 1, iscale = (2**ilevel_max) / 2,

! Number of zones within a sector, assuming fully filled levels, at the finest AMR level.
     1            ntp_global = ntp_zones * ntp_patches * iscale,
     1            nr_global = nr_zones * nr_patches * iscale,
     1            NGhost_global = NGhost * iscale,

! Number of zones within a sector on the base mesh.
     1            ntp_global_base = ntp_zones * ntp_patches,
     1            nr_global_base = nr_zones * nr_patches,

! Max theoretical patches within a sector for "ilevel_max" hierarchy.
! Also, number of patches at base level within a sector and max patches
! actually allowed on a processor.

     1      Max_base_patches_per_sector = ntp_patches**2 * nr_patches,

     1      Max_Theoretical_patches_per_sector =  ! Sum geometric series with ratio 8.
     1      Max_base_patches_per_sector * ( 8**ilevel_max - 1) / 7,

!     1      Max_patches_per_proc = Max_base_patches_per_sector,  ! See "BUILD_BASE_PATCHES".
     1      Max_patches_per_proc = 1,  ! See "BUILD_BASE_PATCHES".

! MHD parameters.
     1            n_cc_components = NVAR + NFLUID,
     1            n_fc_components = 1, n_fc_modes = 1,

#if ( NUMORDER == 2)
     1            n_cc_modes = 4
#elif ( NUMORDER == 3)
     1            n_cc_modes = 10
#elif ( NUMORDER == 4)
     1            n_cc_modes = 20
#endif

     1            )

#if( USE_B0B1_FORMULATION == 1 )
      CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID + 3
     1                                   + SILO_EXTRA_DATA)
#else
      CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID
     1                                   + SILO_EXTRA_DATA)
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Parent-Child Sibling information is stored on each processor so that
! the full composite mesh is easily accessed.

! TODO variable review got to here _________________________________

      INTEGER, BIND(C) ::

! (-1,-1,-1) if the patch is inactive. If active data is present,
! "iactive_FromSector" gives (sector#, patch#, processsor#).

     1   iactive_FromSector ( 6, Max_Theoretical_patches_per_sector, 3),
! TODO SETBY build_base_patches
! TODO USEDBY build_base_patches

! Uses "iactive_FromSector" to build the sector of the neighboring rt,lf,up,dn,tp,bt
! siblings in "isibling_sector_t". Also gives the number of the neighbor
! in "isibling_number_t" and the processor on which this sibling's data
! resides "isibling_proc_t".


     1   isibling_sector_t ( 6, Max_Theoretical_patches_per_sector, 7),
! TODO ANSWERS  Is my neighbor patch on the same sector as I?

     1   isibling_number_t ( 6, Max_Theoretical_patches_per_sector, 7),
! TODO ANSWERS  What is the number of neighbor patch (isector, ipatch, dir)? ....how can I get the number when I have the number...

     1   isibling_proc_t ( 6, Max_Theoretical_patches_per_sector, 7),
! TODO ANSWERS  What is the number of the neighbor patch?

! The desired usage is that we have "Max_patches_per_proc" number of patches
! on each processor. For each such active patch, we want its sector number
! and its logical patch number in that sector.

     1   iactive_FromProc ( Max_patches_per_proc, 2)
! TODO SETBY

! ----------------------------------------------------------------------

! Minimal amounts of coordinate information is stored from which we can
! regenerate the geometry.

! "coordinates_gr_global" contains the coordinates of the vertices of the
! finest mesh that covers the 2D surface of the unit sphere. As a result, it has
! information for all six sectors. This is built once at the start of execution.

! From this information, each patch at each level in the AMR hierarchy can
! extract its own "vert_coords_patch, imask_active_patch" which allows the
! geometry for that patch to be quickly rebuilt.

      REAL, BIND(C) ::  coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

! "rb_save_global, dr_save_global, rc_save_global" gives us radial coordinate
! information on the finest level. Actually, only "rb_save_global" is used to
! build "vert_coords_patch".

      REAL, BIND(C) :: rb_save_global ( 0 - NGhost_global:
     1                      nr_global + NGhost_global)

      REAL, BIND(C), DIMENSION ( 1  - NGhost_global: nr_global + NGhost_global) ::
     1     dr_save_global, rc_save_global

! ----------------------------------------------------------------------

! Minimal amounts of geometry information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor. This is so that all other
! geometry information can be quickly regenerated from it.

! Vertices of each of the patches.
! "vert_coords_amr" for all patches is matched by "vert_coords_patch" for single patches.

      REAL, BIND(C) :: vert_coords_amr
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3, Max_patches_per_proc)

! Mask to decide whether that zone can be touched; some zones are tri-corners.
! "imask_active_amr" for all patches is matched by "imask_active_patch" for single patches.

      INTEGER, BIND(C) :: imask_active_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost, Max_patches_per_proc)

! The stencil info. is only two-dimensional and worth storing.
! "istencil_amr" for all patches is matched by "istencil_patch" for single patches.

      INTEGER, BIND(C) :: istencil_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost: ntp_zones + NGhost, 2, 12,
     1                       Max_patches_per_proc)

! For each patch, keep track whether its four corners are tri-corners or quad-corners.
! "itri_patch_amr" for all patches empties into
! "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU" for the patches.

      INTEGER, BIND(C) :: itri_patch_amr ( 4, Max_patches_per_proc)

! Sector and level information for each patch.
! "isector_ilevel_amr" for all patches is matched by "isector" which is the sector
! number for this patch.

      INTEGER, BIND(C) :: isector_ilevel_amr ( 2, Max_patches_per_proc)

! Which "ix_patch, iy_patch, iz_patch" is this patch?
! "ixyz_patch_amr" for all patches empties into "ix_patch, iy_patch, iz_patch"
! for each patch.

      INTEGER, BIND(C) :: ixyz_patch_amr ( 3, Max_patches_per_proc)

! Is top or bottom patch. i.e., has border at rmin or rmax. Just needed for ADD
! and PAD boundary conditions.

      INTEGER, BIND(C) :: is_bottom_or_top ( 2, Max_patches_per_proc)

! --------------------------------------------------
! The above variables represent the minimal amount of geometry information that
! we need to store in order to rebuild all the geometry variables.
! But we also strategically store a more geometry information that are
! necessary.


! Next, from the above two, all other geometry and quadrature information is built.

      REAL, BIND(C) :: xf_quadpts_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost:  nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: xf_quadwt_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost:  nr_zones + NGhost,
     1   NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: xf_quadnorml_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

      REAL, BIND(C) :: xf_area_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                   1 - NGhost: ntp_zones + NGhost,
     1                   1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc)

      REAL, BIND(C) :: xf_centroid_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ----------

      REAL, BIND(C) :: yf_quadpts_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost:  nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: yf_quadwt_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                     0 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost:  nr_zones + NGhost,
     1   NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: yf_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                        0 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

      REAL, BIND(C) :: yf_area_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                   0 - NGhost: ntp_zones + NGhost,
     1                   1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc)

      REAL, BIND(C) :: yf_centroid_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                       0 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ----------

      REAL, BIND(C) :: zf_quadpts_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost:  nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: zf_quadwt_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost,
     1                     0 - NGhost:  nr_zones + NGhost,
     1   NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: zf_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost: ntp_zones + NGhost,
     1                        0 - NGhost:  nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: zf_area_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                   1 - NGhost: ntp_zones + NGhost,
     1                   0 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc)

      REAL, BIND(C) :: zf_centroid_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                       1 - NGhost: ntp_zones + NGhost,
     1                       0 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ------------------------------
! Integration points for magnetic field in the multid riemann solver.
! These are also the locations where the magnetic field is directly reconstructed.

      REAL, BIND(C) :: xb_quadpts_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost:  nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: xb_quadwt_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost:  nr_zones + NGhost,
     1   NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: xb_quadnorml_amr ( 0 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ----------

      REAL, BIND(C) :: yb_quadpts_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost:  nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: yb_quadwt_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                     0 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost:  nr_zones + NGhost,
     1   NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: yb_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                        0 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ----------

      REAL, BIND(C) :: zb_quadpts_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost:  nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: zb_quadwt_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost,
     1                     0 - NGhost:  nr_zones + NGhost,
     1   NBFLD_QUAD_PTS, Max_patches_per_proc)

      REAL, BIND(C) :: zb_quadnorml_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost: ntp_zones + NGhost,
     1                        0 - NGhost:  nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS, Max_patches_per_proc)

! ------------------------------

      REAL, BIND(C) :: v_scale_length_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc)

! ------------------------------

      REAL, BIND(C) :: v_quadpts_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost:  nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER, Max_patches_per_proc)

      REAL, BIND(C) :: v_quadwt_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                    1 - NGhost: ntp_zones + NGhost,
     1                    1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER, Max_patches_per_proc)

      REAL, BIND(C) :: v_volume_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                    1 - NGhost: ntp_zones + NGhost,
     1                    1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc)

      REAL, BIND(C) :: v_centroid_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost: ntp_zones + NGhost,
     1                      1 - NGhost:  nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! ------------------------------

      REAL, BIND(C) :: Coef_Basis_Moments_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                              1 - NGhost: ntp_zones + NGhost,
     1                              1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes, Max_patches_per_proc)

! ----------------------------------------------------------------------

! Minimal amounts of MHD information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor.

! Each and every "var_name_amr" variable for the fluid variables will have
! a counterpart "var_name" for the patches.

! Zone-centered MHD variables for all patches on this processor.

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: nr_zones  + NGhost,
     1      n_cc_components, Max_patches_per_proc) :: U_amr, U0_amr


#if( USE_B0B1_FORMULATION == 1 )

! the silo read/write variable is sized as "n_cc_components + B0 ( 3) + entropy".
      REAL, BIND(C) ::  U_silo_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1   n_cc_components + 3 + SILO_EXTRA_DATA, Max_patches_per_proc)

#else

! the silo read/write variable is sized as "n_cc_components + entropy".
      REAL, BIND(C) :: U_silo_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1   n_cc_components + SILO_EXTRA_DATA, Max_patches_per_proc)

#endif

      REAL, BIND(C) :: B0field_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, Max_patches_per_proc)

! Zone-centered entropy density variable for all patches on this processor.

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: nr_zones  + NGhost,
     1      Max_patches_per_proc) :: entropygr_amr, entropygr0_amr

! Extra variables for silo output. divergence and entropy, tempfloor flags.

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1   divgr_amr, flag_entprs_amr, flag_Tfloor_amr

! Face-centered magnetic field variables for all patches on this processor.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: nr_zones  + NGhost,
     1                  Max_patches_per_proc) ::
     1   Bxf_amr,  Byf_amr,  Bzf_amr,
     1   Bxf0_amr, Byf0_amr, Bzf0_amr

! Zone-centered flattener for for all patches on this processor.

      REAL, BIND(C) :: flatten_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc)

! Zone-centered Source for for all patches on this processor.

      REAL, BIND(C) :: Source_amr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, Max_patches_per_proc)

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                  Max_patches_per_proc) ::
     1        CGL_fbeta_gr_amr, Pelectron_gr_amr

! ----------

! Next predicted timestep for each patch on this processor.

      REAL, BIND(C) :: dtcurp1_amr ( Max_patches_per_proc)

! Reduced minimum of next predicted timestep from all patches on this processor.

      REAL, BIND(C) :: dtcurp1_reduce

! Min and max of variables for each patch on this processor.

      REAL, BIND(C), DIMENSION ( n_cc_components + 1, Max_patches_per_proc) ::
     1     V_min_amr, V_max_amr

! Reduced min and max variables from all patches on this processor.

      REAL, BIND(C), DIMENSION ( n_cc_components + 1) ::
     1    V_min_reduce, V_max_reduce

! ----------
#if ( ACCURACY_ANALYSIS == 1)

      REAL, BIND(C), DIMENSION ( n_cc_components + 1, Max_patches_per_proc) ::
     1     U_L1_amr, U_Linf_amr, U_zc_Linf_amr

      REAL, BIND(C), DIMENSION ( n_cc_components + 1) ::
     1     U_L1,        U_Linf,        U_zc_Linf,
     1     U_L1_reduce, U_Linf_reduce, U_zc_Linf_reduce

#endif

! ----------------------------------------------------------------------

! Minimal amounts of flux and electric field information that is stored for each of the
! "Max_patches_per_proc" that reside on this processor.
! The defects in these fluxes and electric fields are also stored here.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1      EntFlux_x_amr, EntFlux_y_amr, EntFlux_z_amr,
     1      Defect_EntFlux_x_amr, Defect_EntFlux_y_amr,
     1      Defect_EntFlux_z_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   Max_patches_per_proc) ::
     1      Efx_amr, Efy_amr, Efz_amr,
     1      Defect_Efx_amr, Defect_Efy_amr, Defect_Efz_amr

! Fluctuation terms.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components, Max_patches_per_proc) ::

     1      Fluct_x_Pl_amr, Defect_Fluct_x_Pl_amr,
     1      Fluct_y_Pl_amr, Defect_Fluct_y_Pl_amr,
     1      Fluct_z_Pl_amr, Defect_Fluct_z_Pl_amr,

     1      Fluct_x_Mi_amr, Defect_Fluct_x_Mi_amr,
     1      Fluct_y_Mi_amr, Defect_Fluct_y_Mi_amr,
     1      Fluct_z_Mi_amr, Defect_Fluct_z_Mi_amr

! ----------------------------------------------------------------------
! Buffer arrays that store "0:N_Ghost" layers to the Right and Left of current patch.
! Here we set "n_cc_components + 4" so that we can move entropy and field information
! also in one single message.

      REAL, BIND(C), DIMENSION ( 0: NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_r_amr, U_snd_l_amr

      REAL, BIND(C), DIMENSION ( 0: NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4) :: U_recev_r, U_recev_l

! ----------
! Layers to the Up and Down of current patch.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0: NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_u_amr, U_snd_d_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0: NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components + 4) :: U_recev_u, U_recev_d

! ----------
! Layers to the Top and Bottom of current patch.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost, 0: NGhost,
     1   n_cc_components + 4, Max_patches_per_proc) ::
     1   U_snd_t_amr, U_snd_b_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1  0 - NGhost: ntp_zones + NGhost, 0: NGhost,
     1   n_cc_components + 4) :: U_recev_t, U_recev_b

! ----------------------------------------------------------------------
! Buffer arrays that store one layer of fluxes and electric fields to the
! Right and Left of current patch. Here we set "2*n_cc_components + 3" so that
! we can move entropy flux and "Ey and Ez" in one single message.

! TODO GOT TO HERE__________________________________

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components,
     1                  Max_patches_per_proc) ::
     1   Flux_snd_r_amr,    Flux_snd_l_amr,
     1   Flux_defect_r_amr, Flux_defect_l_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_r_patch, Flux_snd_l_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_recev_r, Flux_recev_l

! ------------------------------
! x-Fluctuations.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: nr_zones + NGhost,
     1                  n_cc_components, Max_patches_per_proc) ::

     1   Fluct_Pl_snd_r_amr,    Fluct_Pl_snd_l_amr,
     1   Fluct_Mi_snd_r_amr,    Fluct_Mi_snd_l_amr,
     1   Fluct_Pl_defect_r_amr, Fluct_Pl_defect_l_amr,
     1   Fluct_Mi_defect_r_amr, Fluct_Mi_defect_l_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_snd_r_patch, Fluct_Pl_snd_l_patch,
     1   Fluct_Mi_snd_r_patch, Fluct_Mi_snd_l_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_recev_r, Fluct_Pl_recev_l,
     1   Fluct_Mi_recev_r, Fluct_Mi_recev_l

! --------------------------------------------------
! Layers to the Up and Down of current patch. "Ez and Ex" saved in that order.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components,
     1                  Max_patches_per_proc) ::
     1   Flux_snd_u_amr,    Flux_snd_d_amr,
     1   Flux_defect_u_amr, Flux_defect_d_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_u_patch, Flux_snd_d_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_recev_u, Flux_recev_d

! ------------------------------
! y-Fluctuations.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: nr_zones + NGhost,
     1                  n_cc_components, Max_patches_per_proc) ::

     1   Fluct_Pl_snd_u_amr,    Fluct_Pl_snd_d_amr,
     1   Fluct_Mi_snd_u_amr,    Fluct_Mi_snd_d_amr,
     1   Fluct_Pl_defect_u_amr, Fluct_Pl_defect_d_amr,
     1   Fluct_Mi_defect_u_amr, Fluct_Mi_defect_d_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_snd_u_patch, Fluct_Pl_snd_d_patch,
     1   Fluct_Mi_snd_u_patch, Fluct_Mi_snd_d_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_recev_u, Fluct_Pl_recev_d,
     1   Fluct_Mi_recev_u, Fluct_Mi_recev_d

! --------------------------------------------------
! Layers to the Top and Bottom of current patch. "Ex and Ey" saved in that order.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  3 + 2*n_cc_components,
     1                  Max_patches_per_proc) ::
     1   Flux_snd_t_amr,    Flux_snd_b_amr,
     1   Flux_defect_t_amr, Flux_defect_b_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_t_patch, Flux_snd_b_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_recev_t, Flux_recev_b

! ------------------------------
! z-Fluctuations.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: nr_zones + NGhost,
     1                  n_cc_components, Max_patches_per_proc) ::

     1   Fluct_Pl_snd_t_amr,    Fluct_Pl_snd_b_amr,
     1   Fluct_Mi_snd_t_amr,    Fluct_Mi_snd_b_amr,
     1   Fluct_Pl_defect_t_amr, Fluct_Pl_defect_b_amr,
     1   Fluct_Mi_defect_t_amr, Fluct_Mi_defect_b_amr

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_snd_t_patch, Fluct_Pl_snd_b_patch,
     1   Fluct_Mi_snd_t_patch, Fluct_Mi_snd_b_patch

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  n_cc_components) ::
     1   Fluct_Pl_recev_t, Fluct_Pl_recev_b,
     1   Fluct_Mi_recev_t, Fluct_Mi_recev_b

! ----------------------------------------------------------------------

! A layer "NGhost" deep of conserved variables (U) and fields will have to
! be messaged in the x-right (XR) direction. Similar for the x-left (XL) direction.
! Similar for the y-up (YU) direction and y-down (YD) direction.
! Similar for the z-top (ZT) direction and z-bottom (ZB) direction.

! We also build the analogous MPI handles for messaging of fluxes.

      INTEGER(KIND=MPI_ADDRESS_KIND), BIND(C) ::
     1        window_size_xy, window_size_z, window_offset

      INTEGER, BIND(C) :: Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
     1        Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
     1        Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle,

     1        Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
     1        Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
     1        Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle,

     1        type_mpi_real  ! Holds the number of bytes for real datatype.

























! ON-PROCESSOR, SINGLE-PATCH GEOMETRY VARIABLES:-

! Local data, regenerated from "coordinates_gr_global" and "rb_save_global", that
! contains the full geometry information that is being processed on this processor.

! First, this data is used to build "vert_coords_patch" and "imask_active_patch".

      REAL, BIND(C) :: vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER, BIND(C) :: imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: nr_zones + NGhost)

! --------------------------------------------------

! Next, from the above two, all other geometry and quadrature information is built.

      REAL, BIND(C) :: xf_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL, BIND(C) :: xf_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL, BIND(C) :: xf_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, BIND(C) :: xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL, BIND(C) :: yf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL, BIND(C) :: yf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL, BIND(C) :: yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, BIND(C) :: yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL, BIND(C) :: zf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL, BIND(C) :: zf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL, BIND(C) :: zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL, BIND(C) :: zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost),
     1     zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

! ------------------------------
! Integration points for magnetic field in the multid riemann solver.
! These are also the locations where the magnetic field is directly reconstructed.

      REAL, BIND(C) :: xb_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL, BIND(C) :: xb_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL, BIND(C) :: xb_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------

      REAL, BIND(C) :: yb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL, BIND(C) :: yb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL, BIND(C) :: yb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------

      REAL, BIND(C) :: zb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL, BIND(C) :: zb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL, BIND(C) :: zb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! ------------------------------

      REAL, BIND(C) :: v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------

      REAL, BIND(C) :: v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

      REAL, BIND(C) :: v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

      REAL, BIND(C) :: v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

      REAL, BIND(C) :: v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL, BIND(C) :: v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

      INTEGER, BIND(C), DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

      REAL, BIND(C) :: Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ------------------------------

      INTEGER, BIND(C) :: istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------

! ON-PROCESSOR, SINGLE-PATCH FLUID VARIABLES, INCLUDING _SKINNY VERSIONS:-

      REAL, BIND(C) :: U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, BIND(C) :: U_skinny ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL, BIND(C) :: entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL, BIND(C) :: entropygr_skinny ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL, BIND(C) :: V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   Bxf_skinny, Byf_skinny, Bzf_skinny

      REAL, BIND(C) :: B0field ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost, n_cc_modes) ::
     1   Ex_zc, Ey_zc, Ez_zc, Wzc

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                                             n_cc_components) ::
     1   Source

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost ) ::
     1        CGL_fbeta_gr, Un9, Pelectron_gr, Un10

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr, workgr

      REAL, BIND(C), DIMENSION ( n_cc_components + 1) :: V_min_local, V_max_local

! Extra variables for silo output; entropy, tempfloor flags.

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost) ::
     1   flag_entprs, flag_Tfloor

! ------------------------------
! Variables used for div-free reconstruction of magnetic fields:-

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           NBFLD_QUAD_PTS) ::
     1   Bxf_quadpts, Byf_quadpts, Bzf_quadpts

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Store extremal speeds from the 1D Riemann solver.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr

! Fluid fluxes, fluctuations and resolved states from the 1D Riemann solver.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                                          n_cc_components) ::
     1      Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1      Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi

      REAL, BIND(C), DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                                 n_cc_components) :: dU_dt

! Entropy flux from the 1D Riemann solver.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

! Electric fields at the edges of the mesh.

      REAL, BIND(C), DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Efx, Efy, Efz







! ---------------
! Subroutines
! ---------------




! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CONTAINS


#include <301riemann_module/310build_base_patches.f>
!
#include <300riemann_module/320load_local_buffers_uxy.f>
#include <300riemann_module/321load_local_buffers_uz.f>
#include <300riemann_module/322load_amr_fluxes_at_bdy.f>
#include <300riemann_module/323load_local_buffers_flux_xy.f>
#include <300riemann_module/324load_local_buffers_flux_z.f>
!
#include <300riemann_module/330fetch_buffers_fill_in_ux.f>
#include <300riemann_module/331fetch_buffers_fill_in_uy.f>
#include <300riemann_module/332fetch_buffers_fill_in_uy_top.f>
#include <300riemann_module/333fetch_buffers_fill_in_uz.f>
#include <300riemann_module/334fetch_buffers_fill_in_flux_x.f>
#include <300riemann_module/335fetch_buffers_fill_in_flux_y.f>
#include <300riemann_module/336fetch_buffers_fill_in_flux_z.f>
!
#include <300riemann_module/340synchronize_efz_quad_corner.f>
#include <300riemann_module/341unpack_defect_fluxes.f>
!
#include <300riemann_module/350wincreate_messaging.f>
#include <300riemann_module/351winfree_messaging.f>


! test code

!include <geomesh_module/test_messaging_error_calculation.f>
!include <geomesh_module/test_xy_messaging.f>
!include <geomesh_module/print_variables_one_patch_per_sector.f>
! fungible: this one calls pad_boundary_r and that depends on Application_CS
!include <geomesh_module/make_global_accuracy_cs.f>
! silo test
!include <geomesh_module/write_silo_one_block.f>



      END MODULE RIEMANN_3D_CS_AMR_SUPPORT

