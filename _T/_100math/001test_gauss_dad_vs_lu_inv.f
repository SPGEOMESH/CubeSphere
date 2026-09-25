



      PROGRAM TEST_GAUSS_DAD_VS_LU_INV

      IMPLICIT NONE

      INTEGER NDim, NMat, NDim_Max, NMat_Max, i, j, k, ifail

      PARAMETER (NDim = 6, NDim_Max = 6, Nmat = 50, NMat_Max = 50)

      REAL timing_info1, timing_info2, time_vec, time_scal, dd,
     1     smallnum, tinynum

      REAL A ( NDim_Max, NDim_Max, NMat_Max),
     1  AInv ( NDim_Max, NDim_Max, NMat_Max),

     1  A_flip ( NMat_Max, NDim_Max, NDim_Max),
     1  AInv_flip ( NMat_Max, NDim_Max, NDim_Max),

     1  A_Ainv ( NDim_Max, NDim_Max, NMat_Max),
     1  Ainv_A ( NDim_Max, NDim_Max, NMat_Max)

      REAL B ( NDim_Max, NDim_Max, NMat_Max),
     1  BInv ( NDim_Max, NDim_Max, NMat_Max),

     1  B_Binv ( NDim_Max, NDim_Max, NMat_Max),
     1  Binv_B ( NDim_Max, NDim_Max, NMat_Max),

     1  B_pt ( NDim_Max, NDim_Max), BInv_pt ( NDim_Max, NDim_Max)

      INTEGER indx_B ( NDim_Max)

! ----------------------------------------------------------------------

      smallnum = 1.0e-30

      tinynum = 1.0e-32

      OPEN ( unit = 6, file = "invtest_gauss.out")

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Testing Gauss_Jordan inv here "
      WRITE ( 6, *)"Using Matrix size = ", NDim , "X", NDim
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of Matrices = ", NMat
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

      call RANDOM_NUMBER(A)

! To make the code cache-friendly, we transpose multidimensional array locations.

      DO i = 1, NMat_Max
      A_flip ( i, 1: NDim_Max, 1: NDim_Max) =
     1    A ( 1: NDim_Max, 1: NDim_Max, i)
      END DO

      B = A

! ----------------------------------------------------------------------
! Make one call to "MatInv_DAD" and time it.

      CALL CPU_TIME ( timing_info1)

      CALL MatInv_DAD ( A_flip, AInv_flip, NDim, NMat,
     1                  NDim_Max, NMat_Max, tinynum)

      CALL CPU_TIME ( timing_info2)

      time_vec = timing_info2 - timing_info1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed the MatInv_DAD subroutine call : ",
     1      time_vec
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! To make the code cache-friendly, we reverse transpose the resulting inverse.

      DO i = 1, NMat_Max
      AInv ( 1: NDim_Max, 1: NDim_Max, i) =
     1    AInv_flip ( i, 1: NDim_Max, 1: NDim_Max)
      END DO

! ----------------------------------------------------------------------
! Make "NMat" number of calls to "ludcmp" and "lubksb" from Numerical Recipes.
! Then time it.

      CALL CPU_TIME ( timing_info1)

      DO k = 1, NMat

      B_pt ( 1: NDim, 1: NDim) = B ( 1: NDim, 1: NDim, k)

        ifail = 0

        CALL ludcmp ( B_pt, NDim, NDim_Max, indx_B, dd, smallnum, ifail)

        BInv_pt = 0.0

        DO i = 1, NDim
        BInv_pt ( i, i) = 1.0
        END DO

        DO j = 1, NDim
          CALL lubksb ( B_pt, NDim, NDim_Max, indx_B, BInv_pt ( 1, j) )
        END DO

        BInv ( 1: NDim, 1: NDim, k) = BInv_pt ( 1: NDim, 1: NDim)

      END DO

      CALL CPU_TIME ( timing_info2)

      time_scal = timing_info2 - timing_info1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed NMat number of calls to ludcmp, lubksb : ",
     1      time_scal
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! ----------------------------------------------------------------------
! Printing A, Ainv alternately

      DO i = 1, NMat

        A_Ainv ( 1:NDim, 1:NDim, i) = MATMUL( A( 1:NDim, 1:NDim, i),
     1                                     Ainv( 1:NDim, 1:NDim, i) )
        Ainv_A ( 1:NDim, 1:NDim, i) = MATMUL( Ainv( 1:NDim, 1:NDim, i),
     1                                     A( 1:NDim, 1:NDim, i) )

        WRITE ( 6, *)"Matrix No. = ", i

        WRITE ( 6, *)"A"
        CALL PRINT_MAT ( NDim, A(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"Ainv"
        CALL PRINT_MAT ( NDim, Ainv(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"A_Ainv"
        CALL PRINT_MAT ( NDim, A_Ainv(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"Ainv_A"
        CALL PRINT_MAT ( NDim, Ainv_A(1:NDim, 1:NDim, i))

        WRITE ( 6, *)"-----------------------------------------------"

      END DO

! ----------------------------------------------------------------------
! Printing B, Binv alternately

      DO i = 1, NMat

        B_Binv ( 1:NDim, 1:NDim, i) = MATMUL( B( 1:NDim, 1:NDim, i),
     1                                   Binv( 1:NDim, 1:NDim, i) )
        Binv_B ( 1:NDim, 1:NDim, i) = MATMUL( Binv( 1:NDim, 1:NDim, i),
     1                                   B( 1:NDim, 1:NDim, i) )

        WRITE ( 6, *)"Matrix No. = ", i

        WRITE ( 6, *)"B"
        CALL PRINT_MAT ( NDim, B(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"Binv"
        CALL PRINT_MAT ( NDim, Binv(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"B_Binv"
        CALL PRINT_MAT ( NDim, B_Binv(1:NDim, 1:NDim, i))
        WRITE ( 6, *)"Binv_B"
        CALL PRINT_MAT ( NDim, Binv_B(1:NDim, 1:NDim, i))

        WRITE ( 6, *)"-----------------------------------------------"

      END DO

! ----------------------------------------------------------------------

      END PROGRAM TEST_GAUSS_DAD_VS_LU_INV



