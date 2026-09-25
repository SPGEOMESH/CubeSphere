


      PROGRAM TEST_GAUSS_MAT_SOL_VS_LUD_SOL

! This program aims to test the speeds of DAD approach to solve a
! set of linear systems given by "MatSolv_DAD" and the "LUD" approach.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NDim, NMat, NDim_Max, NMat_Max, i, j, k, ifail

      PARAMETER (NDim = 6, NDim_Max = 6, Nmat = 50, NMat_Max = 50)

      REAL timing_info1, timing_info2, dd, smallnum,
     1     time_vec, time_scal, tinynum

      REAL A ( NDim_Max, NDim_Max, NMat_Max),

     1  A_flip ( NMat_Max, NDim_Max, NDim_Max)

      REAL R ( NDim_Max, NMat_Max),
     1  Rtmp ( NDim_Max, NMat_Max),
     1  A_X ( NDim_Max, NMat_Max),

     1  R_flip ( NMat_Max, NDim_Max),

     1  BinvR ( NDim_Max, NMat_Max),
     1  B_X ( NDim_Max, NMat_Max)

      REAL B ( NDim_Max, NDim_Max, NMat_Max),

     1  B_pt ( NDim_Max, NDim_Max), BInv_pt ( NDim_Max, NDim_Max),
     1  R_pt ( NDim_Max)

      INTEGER indx_B ( NDim_Max)

! ----------------------------------------------------------------------

      smallnum = 1.0e-30

      tinynum = 1.0e-32

      OPEN ( unit = 6, file = "soltest_gauss.out")


      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Testing Gauss_Jordan Solve here "
      WRITE ( 6, *)"Using Matrix size = ", NDim , "X", NDim
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of Matrices = ", NMat
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! ----------------------------------------------------------------------

! initailize two matrices here

      call RANDOM_NUMBER(A)
      call RANDOM_NUMBER(R)

! To make the code cache-friendly, we transpose multidimensional arrays'
! "A" and "R" locations .

      DO i = 1, NMat_Max
      A_flip ( i, 1: NDim_Max, 1: NDim_Max) =
     1    A ( 1: NDim_Max, 1: NDim_Max, i)

      R_flip ( i, 1: NDim_Max) = R ( 1: NDim_Max, i)
      END DO

      B = A

! ----------------------------------------------------------------------
! Make one call to "MatSolv_DAD" and time it
! when "R_flip" comes out of "MatSolv_DAD" it contains the the flipped solution
! "X" to  AX = R

      CALL CPU_TIME ( timing_info1)

      CALL MatSolv_DAD ( A_flip, R_flip, NDim, NMat,
     1     NDim_Max, NMat_Max, tinynum)

      CALL CPU_TIME ( timing_info2)

      time_vec = timing_info2 - timing_info1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed the MatSolv_DAD  subroutine call : ",
     1      time_vec
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! To make the code cache-friendly, we reverse transpose the resulting inverse.

      DO i = 1, NMat_Max
      Rtmp ( 1: NDim_Max, i) = R_flip ( i, 1: NDim_Max)
      END DO

! ----------------------------------------------------------------------
! Make "NMat" number of calls to "ludcmp" and "lubksb" from Numerical Recipes.
! Then time it.

      CALL CPU_TIME ( timing_info1)

      DO k = 1, NMat

      B_pt ( 1: NDim, 1: NDim) = B ( 1: NDim, 1: NDim, k)

      R_pt ( 1: NDim) = R ( 1: NDim, k)

        ifail = 0

        CALL ludcmp ( B_pt, NDim, NDim_Max, indx_B, dd, smallnum, ifail)

        CALL lubksb ( B_pt, NDim, NDim_Max, indx_B, R_pt )

        BinvR ( 1: NDim, k) = R_pt ( 1: NDim)

      END DO

! The solution is for "BX = R" is obtined as "BinvR"

      CALL CPU_TIME ( timing_info2)

      time_scal = timing_info2 - timing_info1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed NMat number of calls to ludcmp, lubksb : ",
     1      time_scal
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! ----------------------------------------------------------------------
! Printing "AX" and "R" alternately

      DO i = 1, NMat

        A_X ( 1:NDim, i) = MATMUL( A( 1:NDim, 1:NDim, i),
     1                                  Rtmp( 1:NDim, i) )

        WRITE ( 6, *)"Matrix No. = ", i

        WRITE ( 6, *)"A_X"
        CALL PRINT_VECT ( NDim, A_X(1:NDim, i))
        WRITE ( 6, *)"R"
        CALL PRINT_VECT ( NDim, R(1:NDim, i))

        WRITE ( 6, *)"-----------------------------------------------"

      END DO

! ----------------------------------------------------------------------
! Printing "B_X" and "R" alternately

      DO i = 1, NMat

        WRITE ( 6, *)"Matrix No. = ", i

        B_X ( 1:NDim, i) = MATMUL( B( 1:NDim, 1:NDim, i),
     1                                    BinvR( 1:NDim, i) )

        WRITE ( 6, *)"BinvR"
        CALL PRINT_VECT ( NDim, B_X(1:NDim, i))
        WRITE ( 6, *)"R"
        CALL PRINT_VECT ( NDim, R(1:NDim, i))

        WRITE ( 6, *)"-----------------------------------------------"

      END DO

! ----------------------------------------------------------------------

      END PROGRAM TEST_GAUSS_MAT_SOL_VS_LUD_SOL


