




      SUBROUTINE PRINT_VECT( NDim, vector )

      IMPLICIT NONE

      INTEGER NDim

      REAL vector ( NDim)

!      INTEGER i
      REAl(KIND=4) :: vector_4byte( NDim)

      vector_4byte(:) = vector(:)

!      DO i = 1, NDim
        WRITE(6,*) vector_4byte(1:NDim)
!      END DO

      END SUBROUTINE PRINT_VECT




      SUBROUTINE PRINT_MAT( NDim, matrix )

      IMPLICIT NONE

      INTEGER NDim

      REAL matrix ( NDim, NDim)

      INTEGER i
      REAl(KIND=4) :: matrix_4byte( NDim, NDim)

      matrix_4byte( :, :) = matrix( :, :)

      DO i = 1, NDim
        WRITE(6,*) matrix_4byte(i,1:NDim)
      END DO

      END SUBROUTINE PRINT_MAT



