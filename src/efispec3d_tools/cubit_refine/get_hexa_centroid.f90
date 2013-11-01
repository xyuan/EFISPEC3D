subroutine  get_hexa_centroid(xiloc,etloc,zeloc,absci,coloc,xn,yn,zn,x,y,z)

   implicit none

   integer, intent( in) :: xiloc(:)
   integer, intent( in) :: etloc(:)
   integer, intent( in) :: zeloc(:)
   real   , intent( in) :: absci(:,:)
   real   , intent( in) :: coloc(:)
   real   , intent( in) :: xn(:)
   real   , intent( in) :: yn(:)
   real   , intent( in) :: zn(:)
            
   real   , intent(out) :: x
   real   , intent(out) :: y
   real   , intent(out) :: z

   real, parameter      :: ZERO = 0.0
   real                 :: lagrap_geom
   real                 :: local_lagrap1
   real                 :: local_lagrap2
   real                 :: local_lagrap3

   integer              :: i
   integer              :: n_node_geom
   integer              :: n_node_geom_loc

   n_node_geom     = size(xiloc)
   n_node_geom_loc = size(coloc)

   x = 0.0
   y = 0.0
   z = 0.0

   do i = 1,n_node_geom

      local_lagrap1 = lagrap_geom(xiloc(i),ZERO,n_node_geom_loc,coloc,absci)
      local_lagrap2 = lagrap_geom(etloc(i),ZERO,n_node_geom_loc,coloc,absci)
      local_lagrap3 = lagrap_geom(zeloc(i),ZERO,n_node_geom_loc,coloc,absci)

      x = x + local_lagrap1*local_lagrap2*local_lagrap3*xn(i)
      y = y + local_lagrap1*local_lagrap2*local_lagrap3*yn(i)
      z = z + local_lagrap1*local_lagrap2*local_lagrap3*zn(i)

   enddo

   return

end subroutine get_hexa_centroid
!
!
!**********************************************
real function lagrap_geom(i,x,n,c,a)
!**********************************************
!---->compute lagrange polynomial

   implicit none

   integer, intent(in) :: i !local position (i.e. node) in the line [-1:1] where lagrange polynomial = 1 (0 at others nodes)
   real   , intent(in) :: x !abscissa where the polynomial is calculated
   integer, intent(in) :: n !number of nodes in the line [-1:1]
   real   , intent(in) :: c(n)
   real   , intent(in) :: a(n,n)

   integer             :: j

   lagrap_geom = 1.0

   do j = 1,n
      if (j.ne.i) lagrap_geom = lagrap_geom*( (x-c(j))*a(j,i) )
   enddo

end function lagrap_geom
