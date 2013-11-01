subroutine read_topography(fname,x,y,z,dx,dy,x0,y0,nx,ny)

   implicit none

   character(len=90), intent( in) :: fname
   real, allocatable, intent(out) :: x(:)
   real, allocatable, intent(out) :: y(:)
   real, allocatable, intent(out) :: z(:,:)
   real             , intent(out) :: dx
   real             , intent(out) :: dy
   real             , intent(out) :: x0
   real             , intent(out) :: y0
   integer          , intent(out) :: nx
   integer          , intent(out) :: ny

   integer                        :: ix
   integer                        :: iy
   integer                        :: nodata
   character(len=12)              :: ctmp


   open(unit=1,file=trim(adjustl(fname)))

   read(unit=1,fmt=*) ctmp,nx
   read(unit=1,fmt=*) ctmp,ny
   read(unit=1,fmt=*) ctmp,x0
   read(unit=1,fmt=*) ctmp,y0
   read(unit=1,fmt=*) ctmp,dx
   read(unit=1,fmt=*) ctmp,nodata

   dy = dx

   allocate(x(nx),y(ny),z(nx,ny))

   do ix = 1,nx
      x(ix) = real(ix-1)*dx
   enddo

   do iy = 1,ny
      y(iy) = real(ny-1)*dy - real(iy-1)*dy
   enddo

   do iy = 1,ny
      read(unit=1,fmt=*) (z(ix,iy),ix=1,nx)
   enddo

   close(1)

   return

end subroutine read_topography
