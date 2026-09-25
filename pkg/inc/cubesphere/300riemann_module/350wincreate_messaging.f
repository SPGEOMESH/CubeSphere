! contained in RIEMANN_3D_CS_AMR_SUPPORT

















      SUBROUTINE WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
!
! Written By: Sethupathy S.
!
! In this subroutine the data window handles are created for MHD variables
! and the Flux variables in Steps 1 and 2 respectively.

! STEP 1.
! Create windows for message passing for,
!   "Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
!    Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
!    Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle".
! These handles points to,
!   "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr, U_snd_d_amr,"
!   "U_snd_t_amr, U_snd_b_amr".
!
! STEP 2.
! Create windows for message passing for,
!   "Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
!    Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
!    Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle".
! These handles points to,
!   "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr,"
!   "Flux_snd_d_amr, Flux_snd_t_amr, Flux_snd_b_amr".
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: num_of_reals_xy, num_of_reals_z, ierr

! ----------------------------------------------------------------------

! Type of variable being called in bytes is stored in "type_mpi_real".

      CALL MPI_Type_size ( MPI_REAL8, type_mpi_real, ierr)

! ----------------------------------------------------------------------
! STEP 1. Create handles for MHD variables.
! ----------------------------------------------------------------------
! Evaluate number of real variables in the six handles.
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message passing data size is same for "x" and "y" axes
! and it is stored in "num_of_reals_xy". The message passing data size
! along "z" axis is different and it is stored in "num_of_reals_z".

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 ) * Max_patches_per_proc

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 ) * Max_patches_per_proc

! Calculate window size in bytes.

      window_size_xy = type_mpi_real * num_of_reals_xy
      window_size_z  = type_mpi_real * num_of_reals_z

! Create six windows for,
!   "Data_U_XR_Msg_handle, Data_U_XL_Msg_handle,
!    Data_U_YU_Msg_handle, Data_U_YD_Msg_handle,
!    Data_U_ZT_Msg_handle, Data_U_ZB_Msg_handle" to hold values.

! "MPI_INFO_NULL" is something that you just do in MPI.

! The six window handle names will later allow us to get handles on the
! relevant data arrays "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr,
! U_snd_d_amr, U_snd_t_amr, U_snd_b_amr" respectively on the remote arrays.

! It is crucial to release it with "MPI_WIN_FREE".

! x-data
      CALL MPI_WIN_CREATE ( U_snd_r_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_l_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_XL_Msg_handle, ierr)

! y-data
      CALL MPI_WIN_CREATE ( U_snd_u_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_YU_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_d_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_YD_Msg_handle, ierr)

! z-data
      CALL MPI_WIN_CREATE ( U_snd_t_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_ZT_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( U_snd_b_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_U_ZB_Msg_handle, ierr)

! ----------------------------------------------------------------------
! STEP 2. Create handles for Flux variables.
! ----------------------------------------------------------------------
! Evaluate number of real variables in the six handles.
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message passing data size is same for "x" and "y" axes
! and it is stored in "num_of_reals_xy". The message passing data size
! along "z" axis is different and it is stored in "num_of_reals_z".

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 )
     1                * ( 3 + 2*n_cc_components )
     1                * Max_patches_per_proc

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( 3 + 2*n_cc_components)
     1                * Max_patches_per_proc

! Calculate window size in bytes.

      window_size_xy = type_mpi_real * num_of_reals_xy
      window_size_z  = type_mpi_real * num_of_reals_z

! Create six windows for,
!   "Data_Flux_XR_Msg_handle, Data_Flux_XL_Msg_handle,
!    Data_Flux_YU_Msg_handle, Data_Flux_YD_Msg_handle,
!    Data_Flux_ZT_Msg_handle, Data_Flux_ZB_Msg_handle" to hold values.

! "MPI_INFO_NULL" is something that you just do in MPI.

! The six window handle names will later allow us to get handles on the
! relevant data arrays "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr,
! Flux_snd_d_amr, Flux_snd_t_amr, Flux_snd_b_amr" respectively on the
! remote arrays.

! It is crucial to release it with "MPI_WIN_FREE".

! x-data
      CALL MPI_WIN_CREATE ( Flux_snd_r_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_l_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_XL_Msg_handle, ierr)

! y-data
      CALL MPI_WIN_CREATE ( Flux_snd_u_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_d_amr, window_size_xy,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_YD_Msg_handle, ierr)

! z-data
      CALL MPI_WIN_CREATE ( Flux_snd_t_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_ZT_Msg_handle, ierr)

      CALL MPI_WIN_CREATE ( Flux_snd_b_amr, window_size_z,
     1                  type_mpi_real, MPI_INFO_NULL, MPI_COMM_WORLD,
     1                  Data_Flux_ZB_Msg_handle, ierr)

      END SUBROUTINE WINCREATE_MESSAGING


