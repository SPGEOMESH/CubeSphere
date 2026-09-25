



!include <cubesphere/100math/149addmul.f>


      PROGRAM TADD

      IMPLICIT NONE

      REAL x0, y0
      REAL ADD

      x0 = 1.2
      y0 = 3.4

      WRITE(*, *) "Hello", ADD(x0, y0)

#ifdef FOO
      WRITE(*,*) "Foo"
#endif

      STOP
      END


