subroutine compute_z(x_node,y_node,z_node,x_topo,y_topo,z_topo,nx_topo,ny_topo,z_max,z_min)

   implicit none

   real             , intent(in   ) :: x_node
   real             , intent(in   ) :: y_node
   real             , intent(inout) :: z_node
   integer          , intent(in   ) :: nx_topo
   integer          , intent(in   ) :: ny_topo
   real             , intent(in   ) :: x_topo(nx_topo)
   real             , intent(in   ) :: y_topo(ny_topo)
   real             , intent(in   ) :: z_topo(nx_topo,ny_topo)
   real             , intent(in   ) :: z_min
   real             , intent(in   ) :: z_max

   integer                          :: iloc
   integer                          :: jloc
   real                             :: z_topo_max
   real                             :: fac

   call pos(x_topo,nx_topo,x_node,iloc)
   if ( (iloc == 0) .or. (iloc >= nx_topo) ) then
      write(*,*) "error: iloc = ",iloc
      stop
   endif

   call pos(y_topo,ny_topo,y_node,jloc)
   if ( (jloc == 0) .or. (jloc >= ny_topo) ) then
      write(*,*) "error: jloc = ",jloc
      stop
   endif

   call bilinear_interp(x_node,y_node,z_topo_max,x_topo,y_topo,z_topo,nx_topo,ny_topo,iloc,iloc+1,jloc,jloc+1)

   fac = (z_node - z_min)/(z_max - z_min)

   z_node = z_node + z_topo_max * fac

   return

end subroutine compute_z
