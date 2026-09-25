





#if ( GPU_EXECUTION == 0)

      SUBROUTINE MatInv_DAD ( A, AInv, NDim, NMat,
     1                        NDim_Max, NMat_Max, tinynum)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes a dime-a-dozen approach to matrix inversion. The matrices are intended
! to be small. Therefore, Gaussian Elimination with partial pivoting (or partial
! scaled pivoting, if chosen) is used.

! The square matrices all have the same dimension, with the result that
! strong vectorization is applied across matrices. To promote vectorization,
! a cache-friendly arrangement of the matrix elements is used. This arrangement
! gives us the full benefit of cache-prefetching.

! Because all the matrices are drawn from physical problems, it is
! assumed that the pivoting does not have to be very complicated and
! that an inverse always exists.

! There are no guardrails to protect us if the pivot element becomes very
! small or zero; i.e., if one of the matrices becomes singular.

! "NDim" is the dimension of the square matrices.
! "NMat" is the number of square matrices over which we will apply vectorization.

! "NDim_Max" is the dimension of the largest square matrices that we can consider.
! "NMat_Max" is the largest number of such square matrices.

! "tinynum" is about 10^-32

! "A ( NMat_Max, NDim_Max, NDim_Max)" is the input array containing
! "NMat" matrices of dimension "NDimXNDim". It will be destroyed.

! "AInv ( NMat_Max, NDim_Max, NDim_Max)" is the output array. It will eventually
! contain "NMat" matrix inverses of dimension "NDimXNDim".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NDim, NMat, NDim_Max, NMat_Max

      REAL tinynum

      REAL A ( NMat_Max, NDim_Max, NDim_Max),
     1     AInv ( NMat_Max, NDim_Max, NDim_Max)

! ----------------------------------------------------------------------

! Variables past this point are internal to this routine.

      INTEGER i, j, nn, ii

      INTEGER IPiv ( NMat_Max)

      REAL PivNormalize ( NMat_Max, NDim_Max),
     1     ExchangeRow ( NMat_Max, NDim_Max),
     1     RPiv ( NMat_Max)

! ----------------------------------------------------------------------
! STEP 1) Initialize "AInv" to be a set of identity matrices.

      DO i = 1, NDim
      DO j = 1, NDim

!$acc loop vector
        DO ii = 1, NMat
        AInv ( ii, i, j) = 0.0
        END DO

      END DO
      END DO


      DO i = 1, NDim

!$acc loop vector
        DO ii = 1, NMat
        AInv ( ii, i, i) = 1.0
        END DO

      END DO

! ----------------------------------------------------------------------
! Beginning of outer loop which will loop over all the rows.
! Here "nn" is the row of interest; where "nn" will loop from the first to last
! rows of each of the matrices.
! ----------------------------------------------------------------------

      DO nn = 1, NDim  ! Start of outer loop. Working on "nn"th row.

! ----------------------------------------------------------------------
! STEP 2) If partial scaled pivoting is chosen, we have to normalize all rows
! at and after "nn". If we comment out this step we will still get partial
! pivoting, but it will not be partial scaled pivoting.

#if ( 1 == 1)

! ------------------------------
! STEP 2.a) Zero out the normalizing pivot value for "i"th row of "ii"th matrix.

      DO i = nn, NDim  ! Looping over all rows that have not been reduced.

!$acc loop vector
        DO ii = 1, NMat
        PivNormalize ( ii, i) = 0.0
        END DO

      END DO

! ------------------------------
! STEP 2.b) Find maximum absolute pivot value for "i"th row of "ii"th matrix.

      DO i = nn, NDim  ! Looping over all rows that have not been reduced.

        DO j = nn, NDim  ! Looping over all non-zero elements in "i"th row.

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          PivNormalize ( ii, i) = AMAX1 ( PivNormalize ( ii, i),
     1                            ABS ( A ( ii, i, j) ) )

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = nn, NDim".

      END DO  ! End of "DO i = nn, NDim".

! ------------------------------
! STEP 2.c) Find its reciprocal.
! This should be a number that never produces divide by zero.

      DO i = nn, NDim  ! Looping over all rows that have not been reduced.

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.

        PivNormalize ( ii, i) = SIGN ( 1.0, PivNormalize ( ii, i))
     1           / AMAX1 ( tinynum, ABS ( PivNormalize ( ii, i)) )

        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO i = nn, NDim".

! ------------------------------
! STEP 2.d) Use "PivNormalize ( ii, i)" to normalize "i"th row of "ii"th matrix.
! Do this for "A" and "AInv".

! ----------

      DO i = nn, NDim  ! Looping over all rows that have not been reduced.

        DO j = nn, NDim  ! Looping over all non-zero elements in "i"th row.

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          A ( ii, i, j) = PivNormalize ( ii, i) * A ( ii, i, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = nn, NDim".

      END DO  ! End of "DO i = nn, NDim".

! ----------

      DO i = nn, NDim  ! Looping over all rows that have not been reduced.

        DO j = 1, NDim  ! Looping over all elements of AInv in "i"th row.

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          AInv ( ii, i, j) = PivNormalize ( ii, i) * AInv ( ii, i, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = 1, NDim".

      END DO  ! End of "DO i = nn, NDim".

! ----------
#endif

! ----------------------------------------------------------------------
! STEP 3) Get here whether we do scaled partial pivoting or only partial
! pivoting. In this step we find the pivot row and exchange it with the
! "nn" row on which we are currently working.

! ------------------------------
! STEP 3.a) Start with "A ( 1: NMat, nn, nn)" as default choice of pivot element.
! "RPiv" should always be a positive definite number. If it is zero, then that
! matrix has become singular, which is unexpected for a physics-based matrix.

! There are no guardrails in this code if such a thing should happen!

!$acc loop vector
      DO ii = 1, NMat
        IPiv ( ii) = nn
        RPiv ( ii) = ABS ( A ( ii, nn, nn) )
      END DO

! ------------------------------
! STEP 3.b) Find the largest value of "RPiv" from all the rows going from
! "nn + 1" to "NDim" and for all the matrices over the range "1: NMat".
! We search the "nn" column of all these rows for the largest pivot element.

      DO i = nn + 1, NDim  ! Looping over all rows from "nn + 1" onwards.

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.

        IF ( RPiv ( ii) .LT. ABS ( A ( ii, i, nn) ) ) THEN

          RPiv ( ii) = ABS ( A ( ii, i, nn) )
          IPiv ( ii) = i

        END IF

        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO i = nn + 1, NDim".

! ------------------------------
! STEP 3.c) Now exchange "nn"th row for "ii"th matrix with the "IPiv ( ii)" row.
! Do this exchange for just the non-zero parts of "A". But also do this exchange
! for the entire row of "AInv".

! ----------

      DO j = nn, NDim  ! Looping over all non-zero elements in rows of "A" that are exchanged.

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.
          ExchangeRow ( ii, j) = A ( ii, nn, j)
          A ( ii, nn, j) = A ( ii, IPiv ( ii), j)
          A ( ii, IPiv ( ii), j) = ExchangeRow ( ii, j)
        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO j = nn, NDim".

! ----------
! This should be a number that never produces divide by zero.

!$acc loop vector
      DO ii = 1, NMat  ! Looping over all matrices.
        RPiv ( ii) = SIGN ( 1.0, A ( ii, nn, nn) )
     1             / AMAX1 ( tinynum, ABS ( A ( ii, nn, nn) ) )
      END DO  ! End of "DO ii = 1, NMat".

! ----------

      DO j = 1, NDim  ! Looping over all elements in rows of "AInv" that are exchanged.

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.
          ExchangeRow ( ii, j) = AInv ( ii, nn, j)
          AInv ( ii, nn, j) = AInv ( ii, IPiv ( ii), j)
          AInv ( ii, IPiv ( ii), j) = ExchangeRow ( ii, j)
        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO j = 1, NDim".

! ------------------------------
! STEP 3.d) Now divide out the "nn"th rows of "A" and "AInv" by "A ( ii, nn, nn)".
! This reciprocal is already stored in "RPiv ( ii)".
! This ensures that "A ( ii, nn, nn)" becomes unity.

      DO j = nn, NDim  ! Looping over all non-zero elements in "nn"th row of "A".

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.

        A ( ii, nn, j) = RPiv ( ii) * A ( ii, nn, j)

        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO j = nn, NDim".

! ----------

      DO j = 1, NDim  ! Looping over all elements in "nn"th row of "AInv".

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.

        AInv ( ii, nn, j) = RPiv ( ii) * AInv ( ii, nn, j)

        END DO  ! End of "DO ii = 1, NMat".

      END DO  ! End of "DO j = 1, NDim".

! ----------------------------------------------------------------------
! STEP 4) By this point, "A ( ii, nn, nn)" has become unity. Because we want
! an inverse, we zero out all elements "A ( ii, 1: nn-1, nn)" and
! "A ( ii, nn+1: NMat, nn)". This is an elementary row operation applied to "A".
! Exactly the same elementary row operation is applied to "AInv".

      DO i = 1, nn - 1  ! Working on rows that are above the "nn" row.

! ----------

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.
          RPiv ( ii) = A ( ii, i, nn)  ! Save the "(i, nn)" elements from all "A".
        END DO  ! End of "DO ii = 1, NMat".

! ----------

        DO j = nn, NDim  ! Looping over all non-zero elements in "i"th row of "A".

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          A ( ii, i, j) = A ( ii, i, j) - RPiv ( ii) * A ( ii, nn, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = nn, NDim".

! ----------

        DO j = 1, NDim  ! Looping over all elements in "i"th row of "AInv".

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          AInv ( ii, i, j) = AInv ( ii, i, j)
     1                     - RPiv ( ii) * AInv ( ii, nn, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = 1, NDim"

! ----------

      END DO  ! End of "DO i = 1, nn - 1".

! ------------------------------

      DO i = nn + 1, NDim  ! Working on rows that are below the "nn" row.

! ----------

!$acc loop vector
        DO ii = 1, NMat  ! Looping over all matrices.
          RPiv ( ii) = A ( ii, i, nn)  ! Save the "(i, nn)" element from all "A".
        END DO  ! End of "DO ii = 1, NMat".

! ----------

        DO j = nn, NDim  ! Looping over all non-zero elements in "i"th row of "A".

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          A ( ii, i, j) = A ( ii, i, j) - RPiv ( ii) * A ( ii, nn, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = nn, NDim".

! ----------

        DO j = 1, NDim  ! Looping over all elements in "i"th row of "AInv".

!$acc loop vector
          DO ii = 1, NMat  ! Looping over all matrices.

          AInv ( ii, i, j) = AInv ( ii, i, j)
     1                     - RPiv ( ii) * AInv ( ii, nn, j)

          END DO  ! End of "DO ii = 1, NMat".

        END DO  ! End of "DO j = 1, NDim"

! ----------

      END DO  ! End of "DO i = nn + 1, NDim"

! ----------------------------------------------------------------------

      END DO  ! End of outer loop. Working on "nn"th row.

! ----------------------------------------------------------------------
! End of outer loop which will loop over all the rows.
! Here "nn" is the row of interest; where "nn" will loop from the first to last
! rows of each of the matrices.
! ----------------------------------------------------------------------

      END SUBROUTINE MatInv_DAD


#endif
! if ( GPU_EXECUTION == 0)

