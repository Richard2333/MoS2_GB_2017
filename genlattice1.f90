module param

   integer, parameter :: row = 7, crow = 7
   integer, parameter :: sca = row*crow
   integer, parameter :: deg = 3
end module
!include "networklibstandard.f90"
include "libdf1.f90"
program main
   use param
   use biolibs_fixedscale
   implicit none
   integer :: thenet(sca, deg)
   integer :: i, j, k
   real :: x0, y0, x1, y1, nl, ny, zz, nx, m, ab, fab
   open (unit=19, file='latpos')
   thenet = 0
   do i = 0, row - 1
      do j = 2, crow
         call pedge(thenet, sca, deg, i*crow + j, i*crow + j - 1) !:行内连接
      end do
      !call pedge(thenet,sca,deg,i*crow+1,i*crow+crow)
   end do

   do i = 1, row - 1
      do j = 1, crow
         if ((i)*crow + j <= sca) then
            if (mod(i, 2) == 1 .and. mod(j, 2) == 1) then
               call pedge(thenet, sca, deg, i*crow + j, (i - 1)*crow + j)
            end if
            if (mod(i, 2) == 0 .and. mod(j, 2) == 0) then
               call pedge(thenet, sca, deg, i*crow + j, (i - 1)*crow + j)
            end if
         end if
      end do
   end do

   !call deledge(thenet,sca,deg,sca,sca-1)
   !all deledge(thenet,sca,deg,1,2)

   !call disorderzbound(thenet,row,crow,2,0.1)
   open (unit=11, file='latticeb')
   open (unit=13, file='latview1.csv')
   !write the crystal lattice
   do i = 1, sca
!~       do j = 1, deg
         write (11, *) i,(thenet(i, j),j=1,deg) !(thenet(i,j),',',j=1,deg)
!~       end do
   end do
   do i = 1, sca
      call coordi(i, crow, row, x0, y0, fab)
      do j = 1, deg
      if (thenet(i, j) /= 0) then
         call coordi(thenet(i, j), crow, row, x1, y1, fab)
         write (13, *) x0, ',', y0, ',', x1, ',', y1
      end if
      end do

   end do

   do i = 1, sca
      m = 0
      do j = 1, deg
         if (thenet(i, j) /= 0) m = 1
      end do
      if (m /= 0) call coordi(i, crow, row, x0, y0, ab)
      write (19, *) x0, y0, ab
   end do

end program main

subroutine disorderzbound(thenet, row, crow, zl, ds)
   use biolibs_fixedscale
   implicit none
   integer :: thenet(row*crow, 3)
   integer :: row, crow, zl, i, j, k
   real :: ds, s
   call sr1and()
   do i = 0, zl - 1
      do j = 1, crow
         call radm(s)
         if (s <= ds) then
            call remove(thenet, row, crow, i*crow + j)
            write (*, *) i*crow + j
         end if
      end do
   end do

   do i = row - zl, row - 1
      do j = 1, crow
         call radm(s)
         if (s <= ds) then
            call remove(thenet, row, crow, i*crow + j)
            write (*, *) i*crow + j
         end if
      end do
   end do
end subroutine disorderzbound

subroutine remove(thenet, row, crow, n)
   use biolibs_fixedscale
   implicit none
   integer :: thenet(row*crow, 3)
   integer :: row, crow, n1, n2, i, j, maxdeg, n
   maxdeg = 3
   n1 = n
   do i = 1, 3
      n2 = thenet(n, i)
      if (n2 /= 0) then

         call deledge(thenet, row*crow, maxdeg, n1, n2)
      end if
   end do
end subroutine remove

subroutine coordi(i, crow, row, x, y, ab)
   implicit none
   integer :: i, crow, row, nx, ny
   real :: x, y, ab
   nx = i/crow
   ny = mod(i, crow)
   if (ny == 0) then
      nx = nx - 1
      ny = crow
   end if

   x = 1.732*ny/2
   y = 1.5*(nx) + (-1)**(nx + ny + 1)*0.25
   if (mod(nx, 2) == mod(ny, 2)) ab = 0
   if (mod(nx, 2) /= mod(ny, 2)) ab = 1
end subroutine coordi
