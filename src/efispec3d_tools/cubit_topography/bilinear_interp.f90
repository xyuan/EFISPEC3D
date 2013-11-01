subroutine bilinear_interp(xn,yn,z,x,y,array,nx,ny,jx,jxp1,jy,jyp1)

   implicit none

   real   , intent( in)      :: xn
   real   , intent( in)      :: yn
   real   , intent(out)      :: z

   real   , dimension(nx)    :: x
   real   , dimension(ny)    :: y
   real   , dimension(nx,ny) :: array

   integer, intent( in)      :: nx
   integer, intent( in)      :: ny

   integer, intent( in)      :: jx
   integer, intent( in)      :: jxp1
   integer, intent( in)      :: jy
   integer, intent( in)      :: jyp1
                             
   real                      :: s1
   real                      :: s2
   real                      :: s3
   real                      :: s4
   real                      :: zt
   real                      :: zta
   real                      :: ztb
   real                      :: zu
   real                      :: zua
   real                      :: zub
!
!--------------------------------------------------------------------
!
   s1 = array(jx  ,jy  )
   s2 = array(jxp1,jy  )
   s3 = array(jxp1,jyp1)
   s4 = array(jx  ,jyp1)
!
!  find slopes used in interpolation
!  i) x.
   zta = xn - x(jx)
   ztb = x(jxp1) - x(jx)
   zt  = zta/ztb
!
!  ii) y.
   zua = yn - y(jy)
   zub = y(jyp1) - y(jy)
   zu  = zua/zub
!
!  use bilinear interpolation
   z = (1.0-zt)*(1.0-zu)*s1 + zt*(1.0-zu)*s2 + zt*zu*s3 + (1.0-zt)*zu*s4

   return

end subroutine bilinear_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pos(xx,n,x,j)

   implicit none

   integer, intent( in) :: n
   real   , intent( in) :: x
   real   , intent( in) :: xx(n)
   integer, intent(out) :: j

   integer :: jl
   integer :: jm
   integer :: ju
!
!  initialize upper & lower limits.
   jl = 0
   ju = n + 1

   do while ( .true. )

      if ( .not..true. ) then
         return

      elseif ( ju-jl>1 ) then

         jm = (ju+jl)/2
         if ( xx(n)>xx(1) .eqv. x>xx(jm) ) then
            jl = jm
         else
            ju = jm
         endif
!
!        repeat untill the test condition is satisfied.
         cycle
      endif
!
!     set the output.
!FDM
      if (jl >= n) jl = n - 1
      if (jl <= 0) jl = 1
!FDM
      j = jl
      exit

   enddo

end subroutine pos
