! program readdata
! API:
! subroutine INCNAME_8 \in THISFILE
! subroutine INCNAME_SP \in THISFILE

      program readdata

      implicit none

      integer nx, ny
      parameter (nx = nxnum, ny = nynum)

      real*4 dataarr(nx,ny)
      real*8 x(nx), y(ny)
      real*4 r4_x, r4_y
      real*8 xmin, xmax, ymin, ymax
      real*8 xb(0:nx), xc(nx), yb(0:ny), yc(ny), dx, dy
      integer i, j, ll
      real*4 maxdata,mindata

      character*4 arg
      character*8 file1

      call getarg(1, arg)
      file1 = arg//"0001"

      xmin = xmin1
      xmax = xmax1
      ymin = ymin1
      ymax = ymax1

      dx = ( xmax - xmin) / nx
      xb ( 0) = xmin
      dy = ( ymax - ymin) / ny
      yb ( 0) = ymin


      DO i = 1, nx
        xb(i) = xb(i-1)+dx
        xc(i) = 0.5*(xb(i)+xb(i-1))
      END DO

      DO i = 1, ny
        yb(i) = yb(i-1)+dy
        yc(i) = 0.5*(yb(i)+yb(i-1))
      END DO


      open (12, file='filelist', status='new')

      do ll = 1, nslice

      open (10, file=file1, status='old')
      open (11, file=file1//'_rewrite', status='new')

      read (10, *) ((dataarr(i,j), i = 1, nx), j = 1, ny)

      maxdata = dataarr (1,1)
      mindata = dataarr (1,1)
      do i = 1, nx
        do j = 1, ny
          r4_x = xc(i)
          r4_y = yc(j)
          write(11,*) r4_x,r4_y,dataarr(i,j)
          if (dataarr(i,j).gt.maxdata) maxdata = dataarr(i,j)
          if (dataarr(i,j).lt.mindata) mindata = dataarr(i,j)
        enddo
        write(11,*)
      enddo

      close(10)
      close(11)
      write(12,*) file1

      CALL INCNAME_8 ( file1, 1)
      enddo

      if (abs(mindata-maxdata).le.1.e-6) then
        if ( abs(mindata).le.1.e-6) then
          mindata = -1.e-6
          maxdata = 1.e-6
        else
          mindata = mindata - abs ( mindata) * 1.e-6
          maxdata = maxdata + abs ( maxdata) * 1.e-6
        endif
      endif

      open (10, file='min', status='new')
      open (11, file='max', status='new')
      write(10,*) mindata
      write(11,*) maxdata
      close(10)
      close(11)
      close(12)


      end program readdata
