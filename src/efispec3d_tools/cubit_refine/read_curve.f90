subroutine read_curve(f,n,x,y)

   implicit none

   character(len= 90), intent( in) :: f
   integer           , intent(out) :: n
   real, allocatable , intent(out) :: x(:)
   real, allocatable , intent(out) :: y(:)
   real                            :: rtmp

   integer                         :: i
   integer                         :: ios


!
!->first pass to count the number of data
   n = 0
   open(unit=1,file=trim(adjustl(f)))
   do
      read(unit=1,fmt=*,iostat=ios) rtmp
      if (ios /= 0) exit
      n = n + 1
   enddo
   rewind(1)

!
!->memory allocation
   n = n + 1
   allocate(x(n),y(n))

!
!->second pass to read x,y coordinates
   do i = 1,n-1
      read(unit=1,fmt=*) x(i),y(i)
   enddo 
   close(1)

!
!->close the curve
   x(n) = x(1)
   y(n) = y(1)

   return

end subroutine read_curve
