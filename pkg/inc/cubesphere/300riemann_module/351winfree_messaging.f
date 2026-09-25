! contained in RIEMANN_3D_CS_AMR_SUPPORT


      SUBROUTINE WINFREE_MESSAGING ()

! ----------------------------------------------------------------------
!
! Written By: Sethupathy S.
!
! Free up all the MHD handles and Flux handles.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: ierr

! Must deallocate. Otherwise expect memory leaks.

! MHD handles.
      CALL MPI_WIN_FREE ( Data_U_XR_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_XL_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_YU_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_YD_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_ZT_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_U_ZB_Msg_handle, ierr)

! Flux handles.
      CALL MPI_WIN_FREE ( Data_Flux_XR_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_XL_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_YU_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_YD_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_ZT_Msg_handle, ierr)
      CALL MPI_WIN_FREE ( Data_Flux_ZB_Msg_handle, ierr)

      END SUBROUTINE WINFREE_MESSAGING
